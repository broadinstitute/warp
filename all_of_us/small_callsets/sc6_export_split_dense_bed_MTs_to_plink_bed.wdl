version 1.0

# Workflow: sc6_export_split_dense_MTs_to_plink_bed
# Summary: Given a split dense MatrixTable with multi-allelic variants for each bed name,
# create corresponding PLINK bed files by chromosome for each bed file.
## All BED files referenced in this python script are UCSC bed files (as opposed to PLINK bed files).

## Copyright Broad Institute, 2023
##
## This WDL pipeline processes a list of dense MatrixTable with bi-allelic variants for each bed file, to create the corresponding
## PLINK bed files by chromosome for each bed file.
## It's designed to be used with human genomic data in the context of variant analysis and annotation.
##
## Requirements/expectations:
## - **This script CAN ONLY BE RUN in Terra and will NOT work in vanilla cromwell**
## - Human genomic data in Hail MatrixTable format
## - UCSC BED files to specify genomic regions of interest, specifically:
##   - Exome regions
##   - ClinVar regions
##   - ACAF Threshold (Common variants)
## - Parameters for Spark and Hail configuration to optimize performance on cloud platforms
## - The script is tested and expected to run on the hg38 reference genome
##
## This WDL script is optimized for cloud-based execution, with specific parameters for Google Cloud Platform.
##
## This WDL is adopted from an example of running Hail in Cromwell/WDL in Terra
## (https://github.com/broadinstitute/aou-ancestry/tree/main/script/wdl/hail_in_wdl)
## For more details: (https://docs.google.com/document/d/1_OY2rKwZ-qKCDldSZrte4jRIZf4eAw2d7Jd-Asi50KE/edit#heading=h.jhakz81sud23)
##
## LICENSING:
## This script is released under the [Appropriate License] (e.g., MIT, BSD-3, etc.). Users are responsible
## for ensuring they are authorized to run all components of this script. Please consult the relevant
## documentation for licensing details of Hail and other tools used in this pipeline.
##
## For information on tool versions and parameters, refer to the specific Docker containers and Hail
## configuration used in this pipeline.

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow export_split_dense_MTs_to_plink_bed {
    input {
        # Analysis Parameters
        String reference_genome  # Reference genome identifier (e.g., "hg38")
        String version_label  # Label to version the run (e.g., "aou_delta")
        String task_label  # Label to organize output (e.g., "bucket/task_label/output")
        String executor_cores  # Number of cores assigned to each Spark executor
        String executor_memory  # Memory assigned to each Spark executor
        String driver_memory  # Memory assigned to the Spark driver
        String output_gs_url  # Google Cloud Storage path for the output MatrixTable
        File input_split_dense_mts_json_path # File name of the bed dense MatrixTables json file, suffix should be .json
        String skip_beds_names # Bed names to skip if not working on all beds
        String chrs # Chromosomes to work on

        # Cluster Parameters
        # Some of these must be taken from Terra workspace cloud information
        Int num_workers  # Number of workers (per shard) for the Hail cluster
        String gcs_project  # Google project ID for Dataproc
        String gcs_subnetwork_name = 'subnetwork'  # Subnetwork name, set if running in Terra Cromwell
        File submission_script  # Script to run on the cluster
        String region = "us-central1"  # Region for Dataproc, set if running in Terra Cromwell

        # VM Parameters
        String hail_docker = "gcr.io/broad-dsde-methods/aou-auxiliary/hail_dataproc_wdl:0.2.125"  # Docker image with Hail and Google Cloud SDK
    }

    call export_split_dense_MTs_to_plink_bed {
        # Task inputs mirror workflow inputs
        input:
            reference_genome = reference_genome,
            version_label = version_label,
            task_label = task_label,
            executor_cores = executor_cores,
            executor_memory = executor_memory,
            driver_memory = driver_memory,
            gcs_project = gcs_project,
            num_workers = num_workers,
            gcs_subnetwork_name = gcs_subnetwork_name,
            submission_script = submission_script,
            hail_docker = hail_docker,
            region = region,
            output_gs_url = output_gs_url,
            input_split_dense_mts_json_path = input_split_dense_mts_json_path,
            skip_beds_names = skip_beds_names,
            chrs = chrs
                                           }
    output {
        File exome_manifest_bim = export_split_dense_MTs_to_plink_bed.exome_manifest_bim
    }
}

task export_split_dense_MTs_to_plink_bed{
    input {
        # Task-specific inputs with descriptions
        String reference_genome
        String version_label
        String task_label
        String executor_cores
        String executor_memory
        String driver_memory
        String output_gs_url
        File submission_script
        String gcs_project
        String region = "us-central1"
        Int num_workers
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name
        String hail_docker
        File input_split_dense_mts_json_path
        String chrs
        String skip_beds_names
    }

    parameter_meta {
        # Additional metadata for input parameters
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
        set -euxo pipefail

        # Configure Google Cloud account and store the active account name
        gcloud config list account --format "value(core.account)" 1> account.txt

        # Verify Python3 availability in the Docker image
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** Version info: $(python3 -V)"
            # Quick test to ensure Python3 is working
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi

        # Inline Python script to start a Hail Dataproc cluster and submit a job
        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        import re
        from google.cloud import dataproc_v1 as dataproc

        # Function to replace unacceptable characters with '-'
        def sanitize_label(label):
            return re.sub(r'[^a-z0-9-]', '-', label.lower())

        def truncate_cluster_name(name):
            # Ensure the name starts with a lowercase letter and convert to lowercase
            name = name.lower()
            if not name[0].isalpha():
                raise ValueError("Cluster name must start with a letter.")

            # Truncate to a maximum of 10 characters
            truncated = name[:10]

            # Ensure the name does not end with a hyphen
            if truncated[-1] == '-':
                truncated = truncated[:-1] + '0'

            return truncated

        # Sanitize version_label and task_label
        sanitized_version_label = sanitize_label("~{version_label}")
        sanitized_task_label = sanitize_label("~{task_label}")
        sanitized_version_label = truncate_cluster_name(sanitized_version_label)
        sanitized_task_label = truncate_cluster_name(sanitized_task_label)

        # Generate a unique cluster name
        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'{sanitized_version_label}-{sanitized_task_label}-hail-{str(uuid.uuid4())[0:13]}'
        print(f"Cluster Name: {cluster_name}")
        script_path = "~{submission_script}"

        # Read the Google Cloud account name
        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            # Construct the command to start the Hail Dataproc cluster using a multi-line f-string
            cluster_start_cmd = f"""
                hailctl dataproc start --num-workers ~{num_workers}
                --region ~{region} --project ~{gcs_project} --service-account {account}
                --worker-machine-type n1-standard-4 --master-machine-type n1-highmem-16 --max-idle=60m --max-age=1440m
                --subnet=projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}
                --enable-component-gateway
                {cluster_name}
            """

            print("Starting cluster...")

            # Replace newline characters with spaces and remove extra spaces
            cluster_start_cmd = ' '.join(cluster_start_cmd.split())

            print(os.popen(cluster_start_cmd).read())

            # Create a Dataproc client and identify the cluster's staging bucket
            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_staging_bucket = cluster.config.temp_bucket

                    # Submit the Hail job to the Dataproc cluster using a multi-line f-string
                    submit_cmd = f"""
                        gcloud dataproc jobs submit pyspark {script_path}
                        --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account}
                        --driver-log-levels root=WARN --
                        --executor_memory ~{executor_memory} --executor_cores ~{executor_cores}
                        --driver_memory ~{driver_memory} --reference_genome ~{reference_genome}
                        --version_label ~{version_label} --task_label ~{task_label}
                        --input_split_dense_mts_json_path ~{input_split_dense_mts_json_path}
                        --skip_beds_names ~{skip_beds_names} --chrs ~{chrs}
                        --output_gs_url ~{output_gs_url}
                    """

                    print("Running: " + submit_cmd)

                    # Replace newline characters with spaces and remove extra spaces
                    submit_cmd = ' '.join(submit_cmd.split())
                    os.popen(submit_cmd).read()

                    print("Copying results out of staging bucket...")
                    staging_cmd = f'gsutil cp -r ~{output_gs_url}/~{version_label}/~{task_label}/exome/plink_bed/manifest.bim.txt  ~{version_label}.exome.manifest.plink_bed.bim.txt'
                    print(staging_cmd)
                    os.popen(staging_cmd).read()

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            # Clean up: stop the Hail Dataproc cluster
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        File exome_manifest_bim = "~{version_label}.exome.manifest.plink_bed.bim.txt"

    }

    runtime {
        # Runtime settings for the task
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}
