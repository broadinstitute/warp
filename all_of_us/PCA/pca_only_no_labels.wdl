version 1.0

workflow pca_only_no_labels {

    input {
        Array[File] hq_sites_vcf_files
        Array[File] hq_sites_vcf_indices
        String final_output_prefix
        Int num_pcs
        Int? min_vcf_partitions_in
    }

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

    call plot_pca {
        input :
            output_prefix=final_output_prefix,
            pca_tsv=create_hw_pca_training.pca_tsv,
            pc1=1,
            pc2=2
    }

    output {
        File training_pca_labels_ht_tsv = create_hw_pca_training.pca_tsv
        File training_pca_labels_tsv_plots = plot_pca.training_pca_labels_tsv_plot
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
        Int disk_gb = 1000 # 1 TB
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
        disk: "${disk_gb} GB"
        preemptible: num_preemptible_attempts
    }
}

task create_hw_pca_training {
    input {
        File full_bgz
        File full_bgz_index
        String final_output_prefix
        Int num_pcs
        Int? min_vcf_partitions_in
        Int disk_gb = 1000
        Int mem_gb = 240
        Int cpu = 48
        String disk_type = "SSD"
    }

    parameter_meta {
        full_bgz: {localization_optional: true}
        full_bgz_index: {localization_optional: true}
    }

    Int min_vcf_partitions = select_first([min_vcf_partitions_in, 200])

    command <<<
        set -e
        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import hail as hl

        spark_conf={
            # Driver settings
            'spark.driver.memory': '40g',
            'spark.driver.maxResultSize': '20g',

            # Executor settings (12 executors × 4 cores = 48 cores total)
            'spark.executor.memory': '18g',
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

        # Write out the training_pca as a tsv.
        scores_training_export_tsv = scores_training.flatten()
        scores_training_export_tsv.export("~{final_output_prefix}_training_pca.tsv")

        EOF
    >>>

    output {
        File pca_tsv = "~{final_output_prefix}_training_pca.tsv"
    }
    runtime {
        docker: "hailgenetics/hail:0.2.67"
        memory: "${mem_gb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${disk_gb} ${disk_type}"
    }
}

task plot_pca {
    input {
        String output_prefix
        File pca_tsv
        Int pc1
        Int pc2
    }

    command <<<
        set -e

        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import ast

        def check_pc(pc: int) -> None:
            if pc < 1:
                raise ValueError(f'Specified pc was negative or zero: {pc}.  Inputs are 1-indexed.')

        def plot_categorical_points(df, score_col, pc1, pc2, col_category, output_figure_fname):
            pc_x = pc1-1
            pc_y = pc2-1

            # Unique labels and colors
            labels = df[col_category].unique()
            colors = plt.cm.rainbow(np.linspace(0, 1, len(labels)))

            display_labels = {label:f'1k-{label}' for label in labels}

            # Plot each group with a different color
            for label, color in zip(labels, colors):
                subset = df[df[col_category] == label]
                x = [d[pc_x] for d in subset['scores']]
                y = [d[pc_y] for d in subset['scores']]
                plt.scatter(x, y, label=display_labels[label], color=color, s=2)

            plt.title('Categorical Computed PCA')
            plt.xlabel(f'PC{pc1}')
            plt.ylabel(f'PC{pc2}')
            plt.legend()

            plt.savefig(output_figure_fname)

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
        df['scores'] = df['scores'].apply(ast.literal_eval)  # Safer than eval()

        # Add pop labels since plot function expects them
        df['pop_label'] = ["No label"] * len(df)
        df['pop_label']

        # Generate the plot
        plot_categorical_points(df, 'scores', pc1, pc2, 'pop_label', output_figure_basename)

        EOF
    >>>

    output {
        File training_pca_labels_tsv_plot = "~{output_prefix}_~{pc1}_~{pc2}.png"
    }

    runtime {
        docker: "faizanbashir/python-datascience:3.6"
        memory: "16 GB"
        cpu: "2"
        disks: "local-disk 500 HDD"
    }
}