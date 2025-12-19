version 1.0

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

    call ConcatenateChromosomalVcfs {
        input:
            vcf_files=hq_sites_vcf_files,
            vcf_indices=hq_sites_vcf_indices,
            output_vcf_basename=final_output_prefix + "_autosomes.vcf.gz"
    }

    # Train the model on the intersection sites (full version that includes the samples)
    call create_hw_pca_training {
        input:
            full_bgz=ConcatenateChromosomalVcfs.concatenated_vcf,
            full_bgz_index=ConcatenateChromosomalVcfs.concatenated_vcf_idx,
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
        File training_pca_labels_tsv_plots_1_2 = plot_1_2.training_pca_labels_tsv_plot
        File training_pca_labels_tsv_plots_3_4 = plot_3_4.training_pca_labels_tsv_plot
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
        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import hail as hl

        spark_conf = {
            # Driver settings
            'spark.driver.memory': '80g',  # Increased driver memory to 80 GB
            'spark.driver.maxResultSize': '20g',

            # Executor settings (12 executors × 4 cores = 48 cores total)
            'spark.executor.memory': '32g',  # Increased executor memory to 32 GB
            'spark.executor.cores': '4',
            'spark.executor.instances': '12',

            # Parallelism settings
            'spark.sql.shuffle.partitions': '1500',
            'spark.default.parallelism': '1500',

            # Adaptive query execution
            'spark.sql.adaptive.enabled': 'true',
            'spark.sql.adaptive.coalescePartitions.enabled': 'true',

            # Memory management
            'spark.executor.memory.fraction': '0.8',
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
        docker: "python:3.8-slim"
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

        def check_pc(pc: int) -> None:
            if pc < 1:
                raise ValueError(f'Specified pc was negative or zero: {pc}.  Inputs are 1-indexed.')

        def plot_categorical_points(df, score_col, pc1, pc2, col_category, output_figure_fname, pct_var=None):
            # Unique labels and colors
            labels = df[col_category].unique()
            colors = plt.cm.rainbow(np.linspace(0, 1, len(labels)))

            # Plot each group with a different color
            for label, color in zip(labels, colors):
                subset = df[df[col_category] == label]
                x = subset[f'PC{pc1}']  # Updated to use individual PC columns
                y = subset[f'PC{pc2}']  # Updated to use individual PC columns
                plt.scatter(x, y, color=color, s=2, alpha=~{alpha})  # Removed legend and added alpha for transparency

            plt.title('CDRv9 PCA')
            if pct_var is not None:
                pct_x = pct_var.get(pc1)
                pct_y = pct_var.get(pc2)
                xlabel = f"PC{pc1} ({pct_x:.2f}%)" if pct_x is not None else f"PC{pc1}"
                ylabel = f"PC{pc2} ({pct_y:.2f}%)" if pct_y is not None else f"PC{pc2}"
            else:
                xlabel = f"PC{pc1}"
                ylabel = f"PC{pc2}"

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

        # Create output filename
        output_figure_basename = f'{output_prefix}_{pc1}_{pc2}.png'

        # Read and process data
        df = pd.read_csv("~{pca_tsv}", sep="\t")
        pct_variance_df = pd.read_csv("~{pct_variance_tsv}", sep="\t")
        pct_variance_lookup = {
            int(str(label).replace("PC", "")): value
            for label, value in zip(pct_variance_df['PC'], pct_variance_df['Percent_Variance_Explained'])
        }

        # Add pop labels since plot function expects them
        df['pop_label'] = ["No label"] * len(df)

        # Generate the plot
        plot_categorical_points(df, 'scores', pc1, pc2, 'pop_label', output_figure_basename, pct_variance_lookup)

        EOF
    >>>

    output {
        File training_pca_labels_tsv_plot = "~{output_prefix}_~{pc1}_~{pc2}.png"
    }

    runtime {
        docker: "faizanbashir/python-datascience:3.6"
        memory: "${mem_gb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${disk_gb} HDD"
    }
}