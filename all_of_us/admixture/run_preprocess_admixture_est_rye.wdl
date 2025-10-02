version 1.0

# Workflow: run_admixture_est_rye
# Description: This workflow generate input data needed by Rye tool (https://github.com/healthdisparities/rye) to generate admixture estimations.
# These input data are following files (more details in the previous link):
## - eigenvalues: This is computed from genomic data (VCF/BED/PED) using PCA
## - Eigenvectors: This is computed from genomic data (VCF/BED/PED) using PCA
## - pop2group: A population to group mapping file

## Copyright Broad Institute, 2023
##
## This WDL pipeline requires a set of files:
## - Eigenvalues (Hail Table): This is computed by Hail PCA (https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca)
## - Training PCA points (tsv): This is computed by Hail PCA function hwe_normalized_pca (https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca)
##   This is applied on training data like HGDP + 1K data
## - Testing PCA points (tsv): [Specific to AoU data, this will be the ancestry prediction file (https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized#h_01GY7QYC017WWMXFB0YTT233G1)]
##   This is computed by Hail function pc_project (https://hail.is/docs/0.2/experimental/index.html#hail.experimental.pc_project)
##   Projecting testing data on the training data PCA space.
##
## The pipeline, then, outputs the files listed at the top of this file in the formats needed by Rye tool
##
## More details can be found on Rye tool github repo (https://github.com/healthdisparities/rye).
##
##
## LICENSING:
## This script is released under the [Appropriate License] (e.g., MIT, BSD-3, etc.). Users are responsible
## for ensuring they are authorized to run all components of this script. Please consult the relevant
## documentation for licensing details of Hail and other tools used in this pipeline.
##
## For information on tool versions and parameters, refer to the specific Docker containers and
## configuration used in this pipeline.

struct RuntimeAttr {
    Float? mem_gb
    Int? disk_gb
    Int? boot_disk_gb
    Int? max_retries
}

workflow run_preprocess_admixture_est_rye {
    input {
        # Analysis Parameters

        String eigenvalues_url  # Path to the Eigenvalues file
        String training_pca_url  # Path to the training PCA data file
        String ancestry_data_url  # Path to the ancestry predictions file that contains testing data PCA data
        String prefix # A prefix appended to the outputs as a task identifer (e.g. aou_delta)

        # VM Parameters
        Int cpus = 16
        # Docker image with Rye tool installed
        String docker_image = "hailgenetics/hail:0.2.67"
    }
    String pipeline_version= "aou_9.0.0"


    call run_preprocess {
        # Task inputs mirror workflow inputs
        input:
            eigenvalues_url = eigenvalues_url,
            training_pca_url = training_pca_url,
            ancestry_data_url = ancestry_data_url,
            prefix = prefix,
            cpus = cpus,
            docker_image = docker_image
    }
    output {
        File rye_eigenval = run_preprocess.rye_eigenval
        File rye_eigenvec = run_preprocess.rye_eigenvec
        File rye_pop2group = run_preprocess.rye_pop2group
    }
}

task run_preprocess {
    input {
        # Task-specific inputs with descriptions
        String eigenvalues_url
        String training_pca_url
        String ancestry_data_url
        String prefix
        Int cpus
        RuntimeAttr? runtime_attr_override
        String docker_image
    }

    RuntimeAttr runtime_default = object {
                                      # Default runtime attributes
                                      mem_gb: 6.5,
                                      disk_gb: 100,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    command <<<

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import pandas as pd



        eigenvalues_df = pd.read_table(f'~{eigenvalues_url}',header=None)
        training_pca = pd.read_table(f'~{training_pca_url}')
        ancestry_data = pd.read_table(f'~{ancestry_data_url}')

        # Create a new DataFrame for Testing PCA features
        testing_pca = ancestry_data[['research_id', 'pca_features']]

        # Rename the columns to match the training data
        testing_pca.columns = ['s', 'scores']

        def transform_dataframe(input_df):
            # Ensure necessary modules are imported
            import pandas as pd
            from ast import literal_eval

            # Initialize the output DataFrame with specific column names
            output_df = pd.DataFrame()

            # Directly assign 'pop_label' to '#FID' if it exists, else assign an empty string
            output_df['#FID'] = input_df['pop_label'] if 'pop_label' in input_df.columns else ""

            # Assign 's' column from input_df to 'IID' in output_df
            output_df['IID'] = input_df['s']

            # Handle 'scores' column based on its content type
            if isinstance(input_df['scores'].iloc[0], str):
                # Convert string representation of lists to actual lists
                scores_array = input_df['scores'].apply(literal_eval)
            elif isinstance(input_df['scores'].iloc[0], list):
                # If already in list format, use directly
                scores_array = input_df['scores']
            else:
                # Log a warning or raise an error for unexpected types
                raise ValueError("Unknown type for 'scores' column")

            # Create a DataFrame from the scores list
            pc_values = pd.DataFrame(scores_array.tolist(), columns=[f"PC{x}" for x in range(1, len(scores_array.iloc[0]) + 1)])

            # Concatenate the new columns to the output DataFrame
            output_df = pd.concat([output_df, pc_values], axis=1)

            return output_df


        output_training_df = transform_dataframe(training_pca)
        output_testing_df = transform_dataframe(testing_pca)

        merged_df = pd.concat([output_training_df, output_testing_df], ignore_index=True)
        merged_df.to_csv(f'~{prefix}_rye.eigenvec',index=False, sep="\t",na_rep="NA")


        distinct_populations = training_pca["pop_label"].drop_duplicates().reset_index(drop=True)
        distinct_populations_df = pd.DataFrame({"Pop": distinct_populations, "Group": distinct_populations})
        distinct_populations_df = distinct_populations_df[distinct_populations_df["Pop"] != "oth"]

        distinct_populations_df.to_csv(f'~{prefix}_rye.pop2group',index=False, sep="\t")

        eigenvalues_df.to_csv(f'~{prefix}_rye.eigenvalues',index=False, sep="\t",header=False)

        EOF

    >>>


    output {
        File rye_eigenval = "~{prefix}_rye.eigenvalues"
        File rye_eigenvec = "~{prefix}_rye.eigenvec"
        File rye_pop2group = "~{prefix}_rye.pop2group"

    }

    runtime {
        # Runtime settings for the task
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: cpus
        docker: docker_image

    }
}
