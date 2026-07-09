version 1.0
import "../../tasks/wdl/Utilities.wdl" as utils

workflow pca_only_no_labels {

    input {
        Array[File] hq_sites_vcf_files
        Array[File] hq_sites_vcf_indices
        String final_output_prefix
        Int num_pcs
        Int? min_vcf_partitions_in
        Float alpha = 0.18

        # PC pairs to plot. If omitted, defaults to plotting PC1 vs PC2 and PC3 vs PC4.
        # If provided, ONLY the given pairs are plotted (include every pair you want).
        # JSON form: [{"left":5,"right":6},{"left":7,"right":8}]
        Array[Pair[Int, Int]]? pc_pairs

        # Optional ancestry-based subsetting. Provide BOTH together (or neither):
        #   - ancestry_list: TSV (with header) of ancestry predictions covering all samples.
        #     Must contain a 'research_id' column (sample ID) and an 'ancestry_pred_other'
        #     column (ancestry label). Other columns are ignored.
        #   - ancestry: which 'ancestry_pred_other' value to subset to (e.g. "eur")
        File? ancestry_list
        String? ancestry
    }
    String pipeline_version = "beta_0.0.0"

    # check that the number of VCF files matches the number of index files and that both are non-zero. Fail if either is 0 or they are different lengths.
    if (length(hq_sites_vcf_files) != length(hq_sites_vcf_indices) || (length(hq_sites_vcf_files) == 0)) {
        call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
            input:
                message = "Input arrays 'hq_sites_vcf_files' and 'hq_sites_vcf_indices' must be non-empty and have the same length."
        }
    }

    if (length(hq_sites_vcf_files) > 1) {
        call ConcatenateChromosomalVcfs {
            input:
                vcf_files=hq_sites_vcf_files,
                vcf_indices=hq_sites_vcf_indices,
                output_vcf_basename=final_output_prefix + "_autosomes.vcf.gz"
        }
    }
    
    File vcf_bgz = select_first([
        ConcatenateChromosomalVcfs.concatenated_vcf, 
        hq_sites_vcf_files[0] # If only one file, use it directly
    ])
    File vcf_bgz_index = select_first([
        ConcatenateChromosomalVcfs.concatenated_vcf_idx, 
        hq_sites_vcf_indices[0] # If only one file, use it directly
    ])

    # Train the model on the intersection sites (full version that includes the samples)
    call create_hw_pca_training {
        input:
            full_bgz=vcf_bgz,
            full_bgz_index=vcf_bgz_index,
            final_output_prefix=final_output_prefix,
            num_pcs=num_pcs,
            min_vcf_partitions_in=min_vcf_partitions_in,
            ancestry_list=ancestry_list,
            ancestry=ancestry
    }

    call compute_pct_variance {
        input:
            eigenvalues_tsv=create_hw_pca_training.pca_eigenvalues_tsv
    }

    call plot_scree {
        input:
            pct_variance_tsv=compute_pct_variance.pct_variance_explained_tsv,
            output_prefix=final_output_prefix
    }

    # Plot the requested PC pairs (defaulting to 1&2 and 3&4 when none are given).
    Array[Pair[Int, Int]] pairs_to_plot = select_first([pc_pairs, [(1, 2), (3, 4)]])

    scatter (pair in pairs_to_plot) {
        call plot_pca {
            input :
                output_prefix=final_output_prefix,
                pca_tsv=create_hw_pca_training.pca_scores_tsv,
                pc1=pair.left,
                pc2=pair.right,
                pct_variance_tsv=compute_pct_variance.pct_variance_explained_tsv,
                alpha=alpha
        }
    }

    output {
        File training_pca_labels_ht_tsv = create_hw_pca_training.pca_scores_tsv
        File training_pca_eigenvalues_tsv = create_hw_pca_training.pca_eigenvalues_tsv
        File training_pca_scree_plot = plot_scree.training_pca_scree_plot
        # One entry per plotted PC pair (in the order of pc_pairs / the defaults).
        # Individual files remain named by PC pair, e.g. <prefix>_1_2_scatter.png.
        Array[File] training_pca_scatter_plots = plot_pca.training_pca_scatter_plot
        Array[File] training_pca_hexbin_plots = plot_pca.training_pca_hexbin_plot
        Array[File] training_pca_3d_density_interactive_plots = plot_pca.training_pca_3d_density_interactive_plot
    }
}

task ConcatenateChromosomalVcfs {
    input {
        # Array of gzipped VCF files, assumed to be in chromosomal order.
        Array[File] vcf_files

        # Array of corresponding gzipped VCF index files (.tbi).
        Array[File] vcf_indices

        # Desired basename for the output concatenated VCF file.
        String output_vcf_basename = "combined_autosomes.vcf.gz"

        # Runtime parameters
        String bcftools_docker = "mgibio/bcftools-cwl:1.12"
        Int memory_gb = 128
        Int cpu = 16
        Int disk_gb = 1500 # 1.5 TB
        Int num_preemptible_attempts = 1
    }

    command <<<
        # Concatenate the VCF files.
        bcftools concat -Oz -o ~{output_vcf_basename} ~{sep=' ' vcf_files}
        bcftools index -t ~{output_vcf_basename}
    >>>

    output {
        File concatenated_vcf = "${output_vcf_basename}"
        File concatenated_vcf_idx = "${output_vcf_basename}.tbi"
    }

    runtime {
        docker: bcftools_docker
        memory: "${memory_gb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${disk_gb} HDD"
        preemptible: num_preemptible_attempts
        bootDiskSizeGb: 50
    }
}

task create_hw_pca_training {
    input {
        File full_bgz
        File full_bgz_index
        String final_output_prefix
        Int num_pcs
        Int? min_vcf_partitions_in
        File? ancestry_list
        String? ancestry
        Int disk_gb = 2000
        Int mem_gb = 512
        Int cpu = 48
        String disk_type = "SSD"
    }

    parameter_meta {
        full_bgz: {localization_optional: true}
        full_bgz_index: {localization_optional: true}
    }

    Int min_vcf_partitions = select_first([min_vcf_partitions_in, 100])

    command <<<
        set -e

        # Reserve ~50 GB for the OS and Python process; give the rest to the JVM driver heap.
        # spark.driver.memory in spark_conf cannot resize a heap that is already running —
        # SPARK_DRIVER_MEMORY must be set before the JVM starts.
        export SPARK_DRIVER_MEMORY="~{mem_gb - 50}g"

        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import hail as hl

        spark_conf = {
            # Driver heap is controlled by SPARK_DRIVER_MEMORY env var set above.
            # Setting spark.driver.memory here has no effect once the JVM is running.
            'spark.driver.maxResultSize': '50g',

            # Parallelism (local mode uses all cores on the VM)
            'spark.sql.shuffle.partitions': '1500',
            'spark.default.parallelism': '1500',

            # Adaptive query execution
            'spark.sql.adaptive.enabled': 'true',
            'spark.sql.adaptive.coalescePartitions.enabled': 'true',

            # Serialization
            'spark.serializer': 'org.apache.spark.serializer.KryoSerializer'
        }
        hl.init(default_reference='GRCh38', idempotent=False, spark_conf=spark_conf)
        print(hl.spark_context().master)

        def subset_to_ancestry(mt, ancestry_list_path:str, selected_ancestry:str):
            """Subset the MatrixTable columns to the samples of one ancestry.

            Reads the ancestry-prediction table by column name: 'research_id' (sample ID)
            and 'ancestry_pred_other' (ancestry label); all other columns are ignored.
            Hard-fails if the selected ancestry matches no rows in the list, or if any
            listed sample is absent from the MatrixTable (exact sample-name match on 's').
            """
            # Read the ancestry-prediction table and keep only the rows whose
            # 'ancestry_pred_other' value matches the requested ancestry.
            anc = hl.import_table(ancestry_list_path).key_by('research_id')
            anc = anc.filter(anc.ancestry_pred_other == selected_ancestry)
            n_selected = anc.count()
            if n_selected == 0:
                raise ValueError(
                    f"Ancestry '{selected_ancestry}' matched no rows in the ancestry list "
                    f"file '{ancestry_list_path}'. Check the 'ancestry' input and the "
                    f"file's 'ancestry_pred_other' column values."
                )

            # Find listed samples that are absent from the matrix table (exact match on 's').
            mt_samples = mt.cols().key_by('s')
            mt_samples = mt_samples.key_by(research_id=mt_samples.s).select()
            missing = anc.anti_join(mt_samples)
            n_missing = missing.count()
            if n_missing == n_selected:
                raise ValueError(
                    f"None of the {n_selected} sample(s) for ancestry '{selected_ancestry}' "
                    f"were found in the matrix table (zero overlap). Verify that the sample "
                    f"IDs match the VCF sample names exactly."
                )
            if n_missing > 0:
                examples = [r.research_id for r in missing.take(5)]
                raise ValueError(
                    f"{n_missing} of {n_selected} sample(s) for ancestry "
                    f"'{selected_ancestry}' are not present in the matrix table "
                    f"(exact match required). Examples: {examples}"
                )

            # All listed samples are present; keep exactly those columns.
            mt = mt.filter_cols(hl.is_defined(anc[mt.s]))
            print(f"Subset matrix table to {n_selected} sample(s) for ancestry '{selected_ancestry}'.")
            return mt

        def get_PCA_scores(vcf_bgz:str, num_pcs:int, min_vcf_partitions=200, ancestry_list_path="", selected_ancestry=""):
            v = hl.import_vcf(vcf_bgz, force_bgz=True,  min_partitions=min_vcf_partitions)
            # Subset to the requested ancestry's samples BEFORE any training.
            if ancestry_list_path:
                v = subset_to_ancestry(v, ancestry_list_path, selected_ancestry)
            eigenvalues, scores, _ = hl.hwe_normalized_pca(v.GT, k=num_pcs, compute_loadings=False)
            return eigenvalues, scores

        # Optional ancestry subsetting inputs (empty strings when not provided).
        ancestry_list_path = "~{default="" ancestry_list}"
        selected_ancestry = "~{default="" ancestry}"
        if bool(ancestry_list_path) != bool(selected_ancestry):
            raise ValueError("'ancestry_list' and 'ancestry' must be provided together, or neither.")

        eigenvalues_training, scores_training = get_PCA_scores(
            "~{full_bgz}", ~{num_pcs}, ~{min_vcf_partitions}, ancestry_list_path, selected_ancestry)

        # Expand PCA scores into separate columns
        scores_training = scores_training.annotate(**{f'PC{i+1}': scores_training.scores[i] for i in range(~{num_pcs})})
        scores_training_export_tsv = scores_training.key_by().select('s', *[f'PC{i+1}' for i in range(~{num_pcs})])

        # Write out the training_pca as a tsv.
        scores_training_export_tsv.export("~{final_output_prefix}_training_pca.tsv")

        # Write out the eigenvalues as a tsv.
        eigenvalues_df = pd.DataFrame(eigenvalues_training, columns=['eigenvalues'])
        eigenvalues_df.to_csv("~{final_output_prefix}_training_pca_eigenvalues.tsv", sep='\\t', index=False)

        EOF
    >>>

    output {
        File pca_scores_tsv = "~{final_output_prefix}_training_pca.tsv"
        File pca_eigenvalues_tsv = "~{final_output_prefix}_training_pca_eigenvalues.tsv"
    }
    runtime {
        docker: "hailgenetics/hail:0.2.134-py3.11"
        memory: "${mem_gb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${disk_gb} ${disk_type}" # large SSD is recommended for increased processing speed
    }
}

task compute_pct_variance {
    input {
        File eigenvalues_tsv
    }

    command <<<
        set -e
        python3 <<EOF
        import os
        import os.path
        import pandas as pd

        def compute_pct_variance(eigenvalues):
            # NOTE: this normalizes by the sum of the COMPUTED eigenvalues (top
            # num_pcs), not the total variance (trace / all eigenvalues). So each
            # value is a PC's share of the retained eigenvalue mass, which sums to
            # 100% over the selected PCs by construction -- it is NOT the true
            # proportion of total variance explained. Labeled accordingly downstream.
            total = sum(eigenvalues)
            pct_var = [100 * v / total for v in eigenvalues]
            return pct_var

        # Read eigenvalues from TSV
        df = pd.read_csv("~{eigenvalues_tsv}", sep="\\t")
        eigenvalues = df['eigenvalues'].tolist()

        pct_var = compute_pct_variance(eigenvalues)

        # Write out each PC's proportion of variance among the computed PCs (%).
        with open("pct_variance_explained.tsv", "w") as f:
            f.write("PC\tVariance_Proportion_Among_Computed_PCs_Pct\\n")
            for i, pv in enumerate(pct_var):
                f.write(f"PC{i+1}\\t{pv}\\n")

        EOF
    >>>

    output {
        File pct_variance_explained_tsv = "pct_variance_explained.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/warp-tools:2.6.1"
        memory: "16 GB"
        cpu: 2
        disks: "local-disk 250 HDD"
    }
}

task plot_scree {
    input {
        File pct_variance_tsv
        String output_prefix

        # Runtime parameters
        Int disk_gb = 100
        Int mem_gb = 8
        Int cpu = 2
    }

    command <<<
        set -e

        python3 <<EOF
        import pandas as pd
        import matplotlib.pyplot as plt

        output_prefix = "~{output_prefix}"

        # Reuse the metric computed by compute_pct_variance (single source of truth)
        # rather than recomputing from eigenvalues. Values are already the per-PC
        # proportion of variance among the computed PCs, expressed as a percent.
        df = pd.read_csv("~{pct_variance_tsv}", sep="\t")
        pct = df['Variance_Proportion_Among_Computed_PCs_Pct'].tolist()

        # Scree plot vs PC index. The "elbow" where the curve levels off indicates
        # how many PCs capture the meaningful structure.
        plt.figure(figsize=(8, 5))
        plt.plot(
            range(1, len(pct) + 1),
            pct,
            marker='o'
        )
        plt.xlabel("Principal Component")
        plt.ylabel("Proportion of variance among computed PCs (%)")
        plt.title("Scree Plot")
        plt.grid(True)

        # Save a static PNG and close the figure (consistent with the other plots;
        # plt.show() would do nothing in this headless task).
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_scree.png")
        plt.close()

        EOF
    >>>

    output {
        File training_pca_scree_plot = "~{output_prefix}_scree.png"
    }

    runtime {
        docker: "faizanbashir/python-datascience:3.6"
        memory: "${mem_gb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${disk_gb} HDD"
    }
}

task plot_pca {
    input {
        String output_prefix
        File pca_tsv
        Int pc1
        Int pc2
        File pct_variance_tsv
        Float alpha # default is 0.18

        # Runtime parameters
        Int disk_gb = 500
        Int mem_gb = 16
        Int cpu = 2
    }

    command <<<
        set -e

        # plotly is not present in this image; install it for the interactive HTML plot.
        # 5.18.0 is the newest release supporting the image's Python 3.6.
        python3 -m pip install --quiet plotly==5.18.0

        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.ndimage import gaussian_filter
        import plotly.graph_objects as go

        def check_pc(pc: int) -> None:
            if pc < 1:
                raise ValueError(f'Specified pc was negative or zero: {pc}.  Inputs are 1-indexed.')

        def plot_scatter(df, pc1, pc2, col_category, output_figure_fname, pct_var=None):
            """Scatter plot."""
            labels = df[col_category].unique()
            colors = plt.cm.rainbow(np.linspace(0, 1, len(labels)))

            for label, color in zip(labels, colors):
                subset = df[df[col_category] == label]
                x = subset[f'PC{pc1}']
                y = subset[f'PC{pc2}']
                plt.scatter(x, y, color=color, s=2, alpha=~{alpha})

            plt.title('CDRv9 PCA')
            if pct_var is not None:
                pct_x = pct_var.get(pc1)
                pct_y = pct_var.get(pc2)
                xlabel = f"PC{pc1} ({pct_x:.2f}%)" if pct_x is not None else f"PC{pc1}"
                ylabel = f"PC{pc2} ({pct_y:.2f}%)" if pct_y is not None else f"PC{pc2}"
            else:
                xlabel = f"PC{pc1}"
                ylabel = f"PC{pc2}"

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.tight_layout()
            plt.savefig(output_figure_fname)
            plt.close()

        def plot_hexbin(df, pc1, pc2, output_figure_fname, pct_var=None):
            fig, ax = plt.subplots()
            # Use a logarithmic color scale (bins='log' -> color encodes log10(count))
            # so fine structure in low-density regions is visible alongside dense cores.
            # viridis (purple->blue->green->yellow) matches the 3D density plot's colormap.
            hb = ax.hexbin(df[f'PC{pc1}'], df[f'PC{pc2}'], gridsize=100, cmap='viridis', mincnt=1, bins='log')
            plt.colorbar(hb, ax=ax, label='log10(count)')
            plt.title('CDRv9 PCA (Hexbin Density, log scale)')
            if pct_var is not None:
                pct_x = pct_var.get(pc1)
                pct_y = pct_var.get(pc2)
                xlabel = f"PC{pc1} ({pct_x:.2f}%)" if pct_x is not None else f"PC{pc1}"
                ylabel = f"PC{pc2} ({pct_y:.2f}%)" if pct_y is not None else f"PC{pc2}"
            else:
                xlabel = f"PC{pc1}"
                ylabel = f"PC{pc2}"
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_aspect('equal', adjustable='box')
            plt.tight_layout()
            plt.savefig(output_figure_fname)
            plt.close()

        def plot_3d_density_interactive(df, pc1, pc2, output_figure_fname, pct_var=None):
            """Interactive (HTML) 3D surface plot of the smoothed log-density of two PCs.

            Builds a 2D histogram, applies a log-density transform and Gaussian
            smoothing, and renders the surface with plotly so it can be rotated and
            zoomed in a browser.
            """
            # Extract the two PCs of interest from the DataFrame
            x = df[f'PC{pc1}']
            y = df[f'PC{pc2}']

            # Same automatic-but-capped binning as the static version
            max_bins = 200
            nbins_x = min(len(np.histogram_bin_edges(x, bins='auto')) - 1, max_bins)
            nbins_y = min(len(np.histogram_bin_edges(y, bins='auto')) - 1, max_bins)
            counts, xedges, yedges = np.histogram2d(x, y, bins=[nbins_x, nbins_y])

            # Log-density transform: compress dynamic range so dense and sparse
            # regions are both visible on the surface
            Z = np.log1p(counts)

            # Gaussian smoothing (same sigma as the static version)
            Z = gaussian_filter(Z, sigma=1.5)

            # Bin centers for the X and Y axes (edges -> centers)
            xcenters = (xedges[:-1] + xedges[1:]) / 2
            ycenters = (yedges[:-1] + yedges[1:]) / 2

            # Transpose Z so its orientation matches the axes: histogram2d indexes
            # counts as [x, y], while plotly's Surface expects z indexed as z[y][x].
            Z = Z.T

            # Build the interactive surface (Viridis to match the static colormap)
            fig = go.Figure(data=[go.Surface(
                z=Z, x=xcenters, y=ycenters,
                colorscale='Viridis',
                colorbar=dict(title='log(1 + density)'),
            )])

            # Same axis-labeling logic as the static functions
            if pct_var is not None:
                pct_x = pct_var.get(pc1)
                pct_y = pct_var.get(pc2)
                xlabel = f"PC{pc1} ({pct_x:.2f}%)" if pct_x is not None else f"PC{pc1}"
                ylabel = f"PC{pc2} ({pct_y:.2f}%)" if pct_y is not None else f"PC{pc2}"
            else:
                xlabel = f"PC{pc1}"
                ylabel = f"PC{pc2}"

            fig.update_layout(
                title='CDRv9 PCA (3D Density, interactive)',
                scene=dict(
                    xaxis_title=xlabel,
                    yaxis_title=ylabel,
                    zaxis_title='log(1 + density)',
                ),
            )

            # Save as a self-contained interactive HTML file
            fig.write_html(output_figure_fname)

        # Define variables from WDL inputs
        output_prefix = "~{output_prefix}"
        pc1 = ~{pc1}
        pc2 = ~{pc2}

        # Validate PC values
        check_pc(pc1)
        check_pc(pc2)

        # Read and process data
        df = pd.read_csv("~{pca_tsv}", sep="\t")

        # Validate the requested PC pair against the available score columns.
        if pc1 == pc2:
            raise ValueError(f'pc1 and pc2 must differ; both were {pc1}.')
        available_pcs = [c for c in df.columns if c.startswith('PC')]
        for pc in (pc1, pc2):
            if f'PC{pc}' not in df.columns:
                raise ValueError(
                    f'Requested PC{pc} is not present in the scores (available: '
                    f'{", ".join(available_pcs)}). Was num_pcs large enough?'
                )

        pct_variance_df = pd.read_csv("~{pct_variance_tsv}", sep="\t")
        pct_variance_lookup = {
            int(str(label).replace("PC", "")): value
            for label, value in zip(pct_variance_df['PC'], pct_variance_df['Variance_Proportion_Among_Computed_PCs_Pct'])
        }

        # Add pop labels since plot function expects them
        df['pop_label'] = ["No label"] * len(df)

        # 1) Scatter-only plot
        scatter_figure_basename = f'{output_prefix}_{pc1}_{pc2}_scatter.png'
        plot_scatter(df, pc1, pc2, 'pop_label', scatter_figure_basename, pct_variance_lookup)

        # 2) Standalone hexbin density plot
        hexbin_figure_basename = f'{output_prefix}_{pc1}_{pc2}_hexbin.png'
        plot_hexbin(df, pc1, pc2, hexbin_figure_basename, pct_variance_lookup)

        # 3) Interactive 3D density surface plot (HTML)
        density_3d_interactive_basename = f'{output_prefix}_{pc1}_{pc2}_3d_density.html'
        plot_3d_density_interactive(df, pc1, pc2, density_3d_interactive_basename, pct_variance_lookup)

        EOF
    >>>

    output {
        File training_pca_scatter_plot = "~{output_prefix}_~{pc1}_~{pc2}_scatter.png"
        File training_pca_hexbin_plot = "~{output_prefix}_~{pc1}_~{pc2}_hexbin.png"
        File training_pca_3d_density_interactive_plot = "~{output_prefix}_~{pc1}_~{pc2}_3d_density.html"
    }

    runtime {
        docker: "faizanbashir/python-datascience:3.6"
        memory: "${mem_gb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${disk_gb} HDD"
    }
}