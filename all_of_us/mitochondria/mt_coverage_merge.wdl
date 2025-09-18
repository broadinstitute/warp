version 1.0

workflow mt_coverage_merge {

    input {
        # Side inputs with fields that are not incoporated in the samples table
        String coverage_tsv
        String ancestry_tsv
        String dob_tsv


        File full_data_tsv
        File? sample_list_tsv

        ####### DATAPROC CLUSTER PARAMETERS #######
        String gcs_project
        String gcs_subnetwork_name = 'subnetwork'
        String output_bucket_path
        String region = "us-central1"
        String hail_docker = "us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:mtswirl_data_proc_v1"
        String prefix

    }
    String pipeline_version = "1.0.0"

    # CLUSTER PARAMETER INFO #
    parameter_meta {
        gcs_project: "The Google project ID information is necessary when spinning up dataproc. This must match the workspace that this workflow is being run from. eg, 'terra-491d5f31'"
        gcs_subnetwork_name: "Set to 'subnetwork' if running in Terra Cromwell"
        output_bucket_path: "The bucket used to stage input and output files, likely the workspace bucket.  Include 'gs://' prefix and any subdirectory structure, but no trailing '/'"
        region: "Set to 'us-central1' if running in Terra Cromwell"
        hail_docker: "The docker image to be used on the dataproc cluster. This must have both Hail and Google Cloud SDK installed. This also needs MTswirl and its dependencies, which is not included in the default hail docker images."
        prefix: "used to name the data proc cluster. Must be alphanumeric. No special characters."
    }

    if (defined(sample_list_tsv)) {
        call subset_data_table {
            input:
                full_data_tsv = full_data_tsv,
                sample_list_tsv = sample_list_tsv
            }
    }

    File input_table = select_first([subset_data_table.subset_tsv, full_data_tsv])


    call process_tsv_files {
        input:
            coverage_tsv = coverage_tsv,
            ancestry_tsv = ancestry_tsv,
            dob_tsv = dob_tsv,
            input_tsv = input_table
    }

    call annotate_coverage {
        input:
            input_tsv = process_tsv_files.processed_tsv,  # Input TSV file path
            gcs_project = gcs_project,
            gcs_subnetwork_name = gcs_subnetwork_name,
            hail_docker = hail_docker,
            output_bucket_path = output_bucket_path,
            prefix = prefix
    }

    call combine_vcfs {
        input:
            input_tsv = process_tsv_files.processed_tsv,  # Input TSV file path
            coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
            artifact_prone_sites_path = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed",  # Path to artifact-prone sites BED file
            file_name = "combined_vcf.vcf.gz"  # Output file name
    }

    call add_annotations as annotated {
        input:
            coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
            coverage_tsv = process_tsv_files.processed_tsv,  # Path to the coverage input TSV file
            vcf_mt = combine_vcfs.results_tar,  # Path to the MatrixTable
            keep_all_samples = true,
            output_name = "annotated"
    }

    call add_annotations as filt_annotated {
        input:
            coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
            coverage_tsv = process_tsv_files.processed_tsv,  # Path to the coverage input TSV file
            vcf_mt = combine_vcfs.results_tar,  # Path to the MatrixTable
            keep_all_samples = false,
            output_name = "filt_annotated"
    }

    output {
        File processed_tsv = process_tsv_files.processed_tsv
        File output_coverage_ht = annotate_coverage.output_ht
        File combined_vcf = combine_vcfs.results_tar
        File annotated_output_tar = annotated.annotated_output_tar
        File filt_annotated_output_tar = filt_annotated.annotated_output_tar
    }
}

task subset_data_table {
    input {
        File full_data_tsv
        File? sample_list_tsv
    }

    String output_tsv = basename(select_first([sample_list_tsv, full_data_tsv]), ".tsv") + "_data.tsv"

    command <<<
        set -euxo pipefail
        python3 <<'EOF'
        import pandas as pd
        import sys

        # If sample_list_tsv is not defined, just copy the full TSV
        if "~{sample_list_tsv}" == "":
            df_main = pd.read_csv("~{full_data_tsv}", sep="\t", dtype=str)
            df_main.to_csv("~{output_tsv}", sep="\t", index=False)
            sys.exit(0)

        df_main = pd.read_csv("~{full_data_tsv}", sep="\t", dtype=str)
        df_samples = pd.read_csv("~{sample_list_tsv}", sep="\t", header=None, names=["sample_id"], dtype=str)

        # Check if TSV has header (Terra-style: entity:sample_id)
        first_col = df_main.columns[0]
        if first_col.startswith("entity:"):
            id_col = first_col
        else:
            sys.exit("ERROR: Unrecognized format for sample ID column in the full data TSV.")

        df_subset = df_main[df_main[id_col].isin(df_samples["sample_id"])]

        df_subset.to_csv("~{output_tsv}", sep="\t", index=False)
        EOF
    >>>

    output {
        File subset_tsv = "~{output_tsv}"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
        memory: "4 GB"
        cpu: "2"
    }
}

task process_tsv_files {
    input {
        File coverage_tsv  # Path to genomics_metrics_Dec142023_1859_02_tz0000.tsv
        File ancestry_tsv  # Path to echo_v4_r2.ancestry_preds.tsv
        File dob_tsv       # Path to echo_DoB_data.tsv
        File input_tsv     # Input TSV file to process
        String output_tsv_name = "processed_data.tsv"  # Name of the output TSV file
    }

    command <<<
        set -euxo pipefail

        python3 <<EOF

        import pandas as pd
        import numpy as np

        # Load the input TSV into a DataFrame
        df = pd.read_csv("~{input_tsv}", sep="\t")

        # Define only the required columns
        columns_needed = ["research_id", "final_base_level_coverage_metrics", "sample", "final_vcf"]

        # Keep only necessary columns
        filtered_df = df

        # Add and rename columns
        filtered_df["s"] = filtered_df["research_id"]
        filtered_df = filtered_df.rename(columns={
        "research_id": "participant_id",
            "final_base_level_coverage_metrics": "coverage",
            "mt_final_vcf": "final_vcf"
        })
        filtered_df_2 = filtered_df.rename(columns={
            "participant_id": "entity:participant_id",
            "mean_coverage": "mt_mean_coverage",
            "median_coverage": "mt_median_coverage",
            "mt_final_vcf": "final_vcf"
        })
        filtered_df_2["freemix_percentage"] = 0

        # Load additional TSV files
        coveragetsv_df = pd.read_csv("~{coverage_tsv}", sep="\t")
        ancestrytsv_df = pd.read_csv("~{ancestry_tsv}", sep="\t")
        dobtsv_df = pd.read_csv("~{dob_tsv}", sep="\t")

        # Merge with filtered_df on 'research_id' and 's'
        filtered_df_2 = filtered_df_2.merge(
            coveragetsv_df[['research_id', 'mean_coverage', 'biosample_collection_date']],
        left_on='s', right_on='research_id', how='left'
        )
        filtered_df_2.drop(columns=['research_id'], inplace=True)
        filtered_df_2 = filtered_df_2.merge(
            ancestrytsv_df[['research_id', 'ancestry_pred']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df_2.drop(columns=['research_id'], inplace=True)
        filtered_df_2 = filtered_df_2.merge(
            dobtsv_df[['research_id', 'date_of_birth']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df_2.drop(columns=['research_id'], inplace=True)

        # Calculate age
        filtered_df_2['date_of_birth'] = pd.to_datetime(filtered_df_2['date_of_birth'])
        filtered_df_2['biosample_collection_date'] = pd.to_datetime(filtered_df_2['biosample_collection_date'])
        filtered_df_2['age'] = pd.to_numeric(
            np.floor((filtered_df_2['biosample_collection_date'] - filtered_df_2['date_of_birth']).dt.days / 365)
        )
        filtered_df_2['age'] = filtered_df_2['age'].fillna(39).astype(int)

        # Rename columns for compatibility
        filtered_df_2.rename(columns={"mean_coverage": "wgs_mean_coverage"}, inplace=True)
        filtered_df_2.rename(columns={"ancestry_pred": "pop"}, inplace=True)

        # Temporary workaround
        filtered_df_2["wgs_median_coverage"] = filtered_df_2["wgs_mean_coverage"]

        # Filter rows with valid coverage metrics
        filtered_df_2 = filtered_df_2[
            (filtered_df_2['coverage_metrics'].notna()) & (filtered_df_2['coverage_metrics'] != '')
        ]

        # Save the processed DataFrame to a TSV file
        filtered_df_2.to_csv("~{output_tsv_name}", sep="\t", index=False)
        EOF
    >>>

    output {
        File processed_tsv = "~{output_tsv_name}"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
        memory: "4 GB"
        cpu: "2"
    }
}

task annotate_coverage {
    input {
        String input_tsv        # Input TSV file path (optional)
        Int? chunk_size = 100    # Chunk size for combining VCFs (default: 100)
        Boolean overwrite = false  # Overwrite existing files (default: false)
        Boolean keep_targets = false  # Add annotation for target (default: false)
        Boolean hail_only = false  # Skip generating flat files (default: false)
        Int? split_merging = 10  # Number of jobs for splitting merging (default: 1)

        String output_bucket_path # Path to the google bucket for output files, e.g. "gs://my-bucket/my-subdir/". No trailing slash.

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-4"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 2880 # two days
        String gcs_subnetwork_name
        String prefix

        String hail_docker

        # VM runtime attributes - this is for the VM that spins up the dataproc cluster
        Int mem_gb =  8
        Int disk_gb = 16
        Int cpu = 1
        Int preemptible_tries = 0
        Int max_retries =  0
        Int boot_disk_gb = 10
    }

    # define output file urls with file names
    String output_aou_vcf_url = output_bucket_path + "merged_coverage_tsvs.ht"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-hail-step1-{str(uuid.uuid4())[0:13]}'

        # Must be local filepath
        script_path = "/opt/mtSwirl/generate_mtdna_call_mt/Terra/annotate_coverage.py"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                 client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    ~{if overwrite then "--overwrite" else ""} \
                    ~{if keep_targets then "--keep-targets" else ""} \
                    --input-tsv ~{input_tsv} \
                    --output-ht-url ~{output_aou_vcf_url} \
                    --temp-dir gs://{cluster_temp_bucket}/{cluster_name} \
                    --chunk-size ~{chunk_size} \
                    ~{if hail_only then "--hail-only" else ""} \
                    --split-merging ~{split_merging}'''

                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    #command <<<
    #    set -euxo pipefail
#
    #    mkdir -p ./tmp
    #    mkdir -p ./results.ht
#
    #    WORK_DIR=$(pwd)
#
#
    #    # Run the annotate_coverage.py script
    #    python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/annotate_coverage.py \
    #    ~{if overwrite then "--overwrite" else ""} \
    #    ~{if keep_targets then "--keep-targets" else ""} \
    #    --input-tsv ~{input_tsv} \
    #    --output-ht "./merged_coverage_tsvs.ht" \
    #    --temp-dir "./tmp/" \
    #    --chunk-size ~{chunk_size} \
    #    ~{if hail_only then "--hail-only" else ""} \
    #    --split-merging ~{split_merging}
#
#
    #    ## note that both the ht and the mt are outputted by the tool
    #    ls -lh ./merged_coverage_tsvs.ht
    #    ls -lht ./merged_coverage_tsvs.mt
#
    #    # Archive the MatrixTable
    #    tar -czf $WORK_DIR/coverages_tsv.mt.tar.gz ./merged_coverage_tsvs.mt*
    #>>>

    output {
        String output_ht = "coverages_tsv.mt.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0"
        memory: "8 GB"
        cpu: "24"
        disks: "local-disk 100 SSD"
    }
}

task combine_vcfs {
    input {
        String input_tsv        # Input TSV file path
        File coverage_mt_tar    # Tar.gzipped directory of the Hail table
        String a_ref = "GRCh38"           # Reference genome (e.g., GRCh38)
        Boolean overwrite = false  # Overwrite existing files (default: false)
        String vcf_col_name = "final_vcf"     # Column name for VCFs
        String output_bucket = "./results"    # Output bucket path
        String artifact_prone_sites_path = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed"  # Path to artifact-prone sites BED file
        String file_name        # Output file name
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results

        WORK_DIR=$(pwd)

        # Unzip the tar.gz file containing the Hail table
        mkdir -p ./unzipped_coverage.ht
        tar -xzf ~{coverage_mt_tar} -C /cromwell_root/unzipped_coverage.ht
        ls -lh /cromwell_root/unzipped_coverage.ht
        cp -r /cromwell_root/unzipped_coverage.ht/merged_coverage_tsvs.mt /cromwell_root/unzipped_coverage.mt

        # Verify the extracted directory TODO
        if [ ! -d "/cromwell_root/unzipped_coverage.mt" ]; then
            echo "Error: Directory '/cromwell_root/unzipped_coverage.mt' does not exist after extraction."
            exit 1
        fi


        # Run the combine_vcfs.py script
        python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/combine_vcfs.py \
        --input-tsv ~{input_tsv} \
        -c /cromwell_root/unzipped_coverage.mt \
        -a-ref ~{a_ref} \
        --overwrite \
        --vcf-col-name ~{vcf_col_name} \
        --output-bucket ./results \
        --temp-dir ./tmp \
        --artifact-prone-sites-path ~{artifact_prone_sites_path} \
        --file-name ~{file_name} \
        --include-extra-v2-fields


        # Tar zip the results directory
        tar -czf $WORK_DIR/results.tar.gz ./results/combined_vcf.vcf.gz.mt
    >>>

    output {
        File results_tar = "results.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0"
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk 50 SSD"
    }
}

task add_annotations {
    input {
        File coverage_mt_tar    # Tar.gzipped directory of the Hail table
        Boolean keep_all_samples = false  # Keep all samples (default: false)
        String coverage_tsv     # Path to the coverage input TSV file
        File vcf_mt             # Path to the MatrixTable
        String output_name      # directory output name
    }

     command <<<
        set -euxo pipefail

        WORK_DIR=$(pwd)

        # Unzip coverage MatrixTable tarball
        mkdir -p ./unzipped_coverage.mt
        tar -xzf ~{coverage_mt_tar} -C ./unzipped_coverage.mt
        ls -lh ./unzipped_coverage.mt/merged_coverage_tsvs.mt

        # Unzip VCF MatrixTable tarball
        mkdir -p ./unzipped_vcf.mt
        tar -xzf ~{vcf_mt} -C ./unzipped_vcf.mt
        ls -lh ./unzipped_vcf.mt/results/combined_vcf.vcf.gz.mt

        # Verify extraction
        if [ ! -d "./unzipped_coverage.mt" ]; then
            echo "Error: Directory './unzipped_coverage.mt' does not exist after extraction."
            exit 1
        fi

        # Run the add_annotations.py script baked inside mtSwirl clone
        python3 /opt/mtSwirl/generate_mtdna_call_mt/add_annotations.py \
            --sample-stats=~{coverage_tsv} \
            ~{if keep_all_samples then "--keep-all-samples" else ""} \
            --fully-skip-vep \
            --band-aid-dbsnp-path-fix \
            --min-het-threshold 0.05 \
            -v ./~{output_name}/vep \
            -a ~{coverage_tsv} \
            -m unzipped_vcf.mt/results/combined_vcf.vcf.gz.mt \
            -d ./~{output_name}

        # Compress the annotated output directory
        tar -czf $WORK_DIR/annotated_output.tar.gz ~{output_name}
    >>>

    output {
        File annotated_output_tar = "annotated_output.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0"
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk 50 SSD"
    }
}