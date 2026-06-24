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
            min_vcf_partitions_in=min_vcf_partitions_in
    }

    call compute_pct_variance {
        input:
            eigenvalues_tsv=create_hw_pca_training.pca_eigenvalues_tsv
    }

    call plot_pca as plot_1_2 {
        input :
            output_prefix=final_output_prefix,
            pca_tsv=create_hw_pca_training.pca_scores_tsv,
            pc1=1,
            pc2=2, 
            pct_variance_tsv=compute_pct_variance.pct_variance_explained_tsv, 
            alpha=alpha
    }

    call plot_pca as plot_3_4 {
        input :
            output_prefix=final_output_prefix,
            pca_tsv=create_hw_pca_training.pca_scores_tsv,
            pc1=3,
            pc2=4, 
            pct_variance_tsv=compute_pct_variance.pct_variance_explained_tsv, 
            alpha=alpha
    }

    output {
        File training_pca_labels_ht_tsv = create_hw_pca_training.pca_scores_tsv
        File training_pca_eigenvalues_tsv = create_hw_pca_training.pca_eigenvalues_tsv
        File training_pca_scatter_plots_1_2 = plot_1_2.training_pca_scatter_plot
        File training_pca_scatter_plots_3_4 = plot_3_4.training_pca_scatter_plot
        File training_pca_hexbin_plots_1_2 = plot_1_2.training_pca_hexbin_plot
        File training_pca_hexbin_plots_3_4 = plot_3_4.training_pca_hexbin_plot
        File training_pca_3d_density_plots_1_2 = plot_1_2.training_pca_3d_density_plot
        File training_pca_3d_density_plots_3_4 = plot_3_4.training_pca_3d_density_plot
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

        def get_PCA_scores(vcf_bgz:str, num_pcs:int, min_vcf_partitions=200):
            v = hl.import_vcf(vcf_bgz, force_bgz=True,  min_partitions=min_vcf_partitions)
            eigenvalues, scores, _ = hl.hwe_normalized_pca(v.GT, k=num_pcs, compute_loadings=False)
            return eigenvalues, scores

        eigenvalues_training, scores_training = get_PCA_scores("~{full_bgz}", ~{num_pcs}, ~{min_vcf_partitions})

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
            total = sum(eigenvalues)
            pct_var = [100 * v / total for v in eigenvalues]
            return pct_var

        # Read eigenvalues from TSV
        df = pd.read_csv("~{eigenvalues_tsv}", sep="\\t")
        eigenvalues = df['eigenvalues'].tolist()

        pct_var = compute_pct_variance(eigenvalues)

        # Write out percent variance explained
        with open("pct_variance_explained.tsv", "w") as f:
            f.write("PC\tPercent_Variance_Explained\\n")
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

        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.ndimage import gaussian_filter

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
            hb = ax.hexbin(df[f'PC{pc1}'], df[f'PC{pc2}'], gridsize=100, cmap='Blues', mincnt=1)
            plt.colorbar(hb, ax=ax, label='Count')
            plt.title('CDRv9 PCA (Hexbin Density)')
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

        def plot_3d_density(df, pc1, pc2, output_figure_fname, pct_var=None):
            """3D surface plot of the smoothed log-density of two PCs."""
            # Extract the two PCs of interest from the DataFrame
            x = df[f'PC{pc1}']
            y = df[f'PC{pc2}']

            # Choose bin counts automatically with numpy's Freedman-Diaconis rule
            # ('auto') rather than hardcoding. On 500k+ samples that rule can return
            # thousands of bins, producing an unwieldy surface mesh, so we cap each
            # axis at a sane maximum for a readable 3D surface.
            max_bins = 200
            nbins_x = min(len(np.histogram_bin_edges(x, bins='auto')) - 1, max_bins)
            nbins_y = min(len(np.histogram_bin_edges(y, bins='auto')) - 1, max_bins)
            counts, xedges, yedges = np.histogram2d(x, y, bins=[nbins_x, nbins_y])

            # Log-density transform: compresses the dynamic range so dense and
            # sparse regions are both visible on the surface
            Z = np.log1p(counts)

            # Smooth the log-density so the surface is continuous rather than jagged
            Z = gaussian_filter(Z, sigma=1.5)

            # Bin centers for the X and Y axes (edges -> centers)
            xcenters = (xedges[:-1] + xedges[1:]) / 2
            ycenters = (yedges[:-1] + yedges[1:]) / 2

            # Build a meshgrid; transpose Z so its orientation matches X/Y
            # (histogram2d indexes Z as [x, y], surface plots expect [y, x])
            X, Y = np.meshgrid(xcenters, ycenters)
            Z = Z.T

            # Create the 3D surface plot
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0, antialiased=True)
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='log(1 + density)')

            ax.set_title('CDRv9 PCA (3D Density)')

            # Same axis-labeling logic as the other plotting functions
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
            ax.set_zlabel('log(1 + density)')

            # Save a static PNG and close the figure (consistent with the rest of the script)
            plt.tight_layout()
            plt.savefig(output_figure_fname)
            plt.close()

        # Define variables from WDL inputs
        output_prefix = "~{output_prefix}"
        pc1 = ~{pc1}
        pc2 = ~{pc2}

        # Validate PC values
        check_pc(pc1)
        check_pc(pc2)

        # Read and process data
        df = pd.read_csv("~{pca_tsv}", sep="\t")
        pct_variance_df = pd.read_csv("~{pct_variance_tsv}", sep="\t")
        pct_variance_lookup = {
            int(str(label).replace("PC", "")): value
            for label, value in zip(pct_variance_df['PC'], pct_variance_df['Percent_Variance_Explained'])
        }

        # Add pop labels since plot function expects them
        df['pop_label'] = ["No label"] * len(df)

        # 1) Scatter-only plot
        scatter_figure_basename = f'{output_prefix}_{pc1}_{pc2}_scatter.png'
        plot_scatter(df, pc1, pc2, 'pop_label', scatter_figure_basename, pct_variance_lookup)

        # 2) Standalone hexbin density plot
        hexbin_figure_basename = f'{output_prefix}_{pc1}_{pc2}_hexbin.png'
        plot_hexbin(df, pc1, pc2, hexbin_figure_basename, pct_variance_lookup)

        # 3) 3D smoothed log-density surface plot
        density_3d_figure_basename = f'{output_prefix}_{pc1}_{pc2}_3d_density.png'
        plot_3d_density(df, pc1, pc2, density_3d_figure_basename, pct_variance_lookup)

        EOF
    >>>

    output {
        File training_pca_scatter_plot = "~{output_prefix}_~{pc1}_~{pc2}_scatter.png"
        File training_pca_hexbin_plot = "~{output_prefix}_~{pc1}_~{pc2}_hexbin.png"
        File training_pca_3d_density_plot = "~{output_prefix}_~{pc1}_~{pc2}_3d_density.png"
    }

    runtime {
        docker: "faizanbashir/python-datascience:3.6"
        memory: "${mem_gb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${disk_gb} HDD"
    }
}