version 1.0

workflow PCA {

    input {
        File    full_bgz
        File    full_bgz_index
        File    hgdp_metadata_file
        String  final_output_prefix
        Int     num_pcs
        Int?    min_vcf_partitions_in
    }

    # Train the model on the intersection sites (full version that includes the samples)
    call create_hw_pca_training {
        input:
            full_bgz                =   full_bgz,
            full_bgz_index          =   full_bgz,
            hgdp_metadata_file      =   hgdp_metadata_file,
            final_output_prefix     =   final_output_prefix,
            num_pcs                 =   num_pcs,
            min_vcf_partitions_in   =   min_vcf_partitions_in
    }

    call plot_pca {
        input :
            output_prefix               =   final_output_prefix,
            training_pca_labels_ht_tsv  =   create_hw_pca_training.training_pca_labels_ht_tsv
    }

    output {
        File training_pca_labels_ht_tsv     =   create_hw_pca_training.training_pca_labels_ht_tsv
        File training_pca_labels_tsv_plot   =   plot_pca.training_pca_labels_tsv_plot
    }
}

task create_hw_pca_training {
    input {
        File    full_bgz
        File    full_bgz_index
        String  final_output_prefix
        File    hgdp_metadata_file
        Int     num_pcs
        Int?    min_vcf_partitions_in
    }

    parameter_meta {
        full_bgz:       {localization_optional: true}
        full_bgz_index: {localization_optional: true}
    }

    Int min_vcf_partitions  =   select_first([min_vcf_partitions_in, 200])

    # You do not have to use heredoc.  You can also send a python script as a parameter and run it in the command block.
    command <<<
        set -e
        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import hail as hl

        hl.init(default_reference='GRCh38', idempotent=True)

        # Read in metadata for ALL samples
        metadata_pd = pd.read_csv("~{hgdp_metadata_file}", sep="\t")

        def get_PCA_scores(vcf_bgz:str, min_vcf_partitions):
            v = hl.import_vcf(vcf_bgz, force_bgz=True, min_partitions=min_vcf_partitions)
            eigenvalues, scores, loadings = hl.hwe_normalized_pca(v.GT, k=~{num_pcs}, compute_loadings=True)
            return eigenvalues, scores, loadings

        def collapse_fin_to_eur(pop):
            if pop=='fin' or pop=='nfe':
                return 'eur'
            else:
                return pop

        eigenvalues_training, scores_training, loadings_training = get_PCA_scores("~{full_bgz}", ~{min_vcf_partitions})

        # Apply any custom processing to the population labels from the training data
        pop_label_pd = metadata_pd[['s', 'population_inference.pop']]
        pop_label_pd['pop_label'] = metadata_pd['population_inference.pop'].apply(collapse_fin_to_eur)

        # Join the labels to the training PCA feature set.
        pop_label_ht = hl.Table.from_pandas(pop_label_pd).key_by('s')
        training_pca_ht = scores_training.join(pop_label_ht, how='inner')

        # Write out the training_pca as a tsv.
        training_pca_ht_export_tsv = training_pca_ht.flatten()
        training_pca_ht_export_tsv.export("~{final_output_prefix}_training_pca.tsv")

        EOF
    >>>

    runtime {
        docker: "hailgenetics/hail:0.2.67"
        memory: "123 GB"
        cpu: "16"
        disk: "local-disk 500 HDD"
    }

    output {
        File training_pca_labels_ht_tsv =   "~{final_output_prefix}_training_pca.tsv"
    }
}

task plot_pca {
    input {
        String  output_prefix
        File    training_pca_labels_ht_tsv
    }

    command <<<
        set -e
        python3 <<EOF
        import os
        import os.path
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt

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

        pc1 = 1
        pc2 = 2
        output_figure_basename = f'~{output_prefix}_{pc1}_{pc2}.png'
        df = pd.read_csv("~{training_pca_labels_ht_tsv}", sep="\t")
        df['scores'] = df['scores'].apply(eval)
        plot_categorical_points(df, 'scores', pc1, pc2, 'pop_label', output_figure_basename)

        EOF
    >>>

    runtime {
        docker: "faizanbashir/python-datascience:3.6"
        memory: "16 GB"
        cpu: "2"
        disk: "local-disk 500 HDD"
    }

    output {
        File training_pca_labels_tsv_plot   =   "~{output_prefix}_1_2.png"
    }
}