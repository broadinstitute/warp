version 1.0

# Workflow: run_relatedness
# Description: This workflow generates two lists:
# 1- A list of all pairs of samples w/ a kinship score greater than a certain threshold (default: 0.1)
# 2- A list of samples that would need to be removed to remove relatedness cofounds from the full cohort

## Copyright Broad Institute, 2023
##
## This WDL pipeline processes a set of files:
## - HQ_variants.vcf: Full VCF with genotypes of the full cohort. This is produced by AoU ancestry workflow (https://github.com/broadinstitute/aou-ancestry)
## - PCA scores of HQ_variants.vcf: This is produced by AoU ancestry workflow (https://github.com/broadinstitute/aou-ancestry)
##
## The pipeline, then, generates two lists: one is the list related samples within the full cohort.
## The other is a list of samples that can be removed from the full cohort to remove relatedness cofounds
##
## Requirements/expectations:
## - A vcf file of high-quality variants for the full cohort
## - PCA scores based on this vcf file
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
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}


workflow run_relatedness {
    input {
        # Analysis Parameters
        String vcf_url  # Path to VCF file with genotypes
        String pca_scores_url  # Path to PCA scores of the samples within vcf_url
        String task_identifier  # An identifier associated with the run like "aou_delta"
        String statistics = 'kin' # Set of statistics to compute. Default to 'kin' which is the kinship statistic
        Float min_individual_maf = 0.01 #  The minimum individual-specific minor allele frequency.
        Int block_size = 512 # Block size of block matrices used in the algorithm
        Float min_kinship = 0.1 # Pairs of samples with kinship lower than min_kinship are excluded from the results.
        Int min_partitions = 20000  # Minimum number of partitions for optimization
        String gcs_output_url  # Google Cloud Storage path for the output files
        String executor_cores  # Number of cores assigned to each Spark executor
        String driver_cores  # Number of cores assigned to the Spark driver
        String executor_memory  # Memory assigned to each Spark executor
        String driver_memory  # Memory assigned to the Spark driver
        String reference_genome  # Reference genome identifier (e.g., "hg38")


        # Cluster Parameters
        # Some of these must be taken from Terra workspace cloud information
        Int num_workers  # Number of workers (per shard) for the Hail cluster
        String gcs_project  # Google project ID for Dataproc
        String gcs_subnetwork_name = 'subnetwork'  # Subnetwork name, set if running in Terra Cromwell
        File submission_script  # Script to run on the cluster
        String region = "us-central1"  # Region for Dataproc, set if running in Terra Cromwell

        # VM Parameters
        String hail_docker = "us.gcr.io/broad-dsde-methods/lichtens/hail_dataproc_wdl:1.1"  # Docker image with Hail and Google Cloud SDK
        #String hail_docker = "gcr.io/broad-dsde-methods/aou-auxiliary/hail_dataproc_wdl:0.2.125"  # Docker image with Hail and Google Cloud SDK
    }
    String pipeline_version="aou_9.1.0"

    call run_relatedness_task {
        # Task inputs mirror workflow inputs
        input:
            vcf_url = vcf_url,
            pca_scores_url = pca_scores_url,
            task_identifier=task_identifier,
            statistics = statistics,
            min_individual_maf = min_individual_maf,
            block_size = block_size,
            min_kinship = min_kinship,
            gcs_output_url = gcs_output_url,
            min_partitions = min_partitions,
            executor_cores = executor_cores,
            driver_cores =  driver_cores,
            executor_memory = executor_memory,
            driver_memory = driver_memory,
            reference_genome=reference_genome,
            gcs_project = gcs_project,
            num_workers = num_workers,
            gcs_subnetwork_name = gcs_subnetwork_name,
            submission_script = submission_script,
            hail_docker = hail_docker,
            region = region,
    }

    output {
        File relatedness = run_relatedness_task.relatedness
        File relatedness_flagged_samples = run_relatedness_task.relatedness_flagged_samples
    }
}

task run_relatedness_task {
    input {
        # Task-specific inputs with descriptions
        String vcf_url
        String pca_scores_url
        String task_identifier
        String statistics
        Float min_individual_maf
        Int block_size
        Float min_kinship
        Int min_partitions
        String gcs_output_url
        String executor_cores
        String executor_memory
        String driver_cores
        String driver_memory
        String reference_genome
        File submission_script
        String gcs_project
        String region
        Int num_workers
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name
        String hail_docker
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

        # Sanitize task_identifier
        sanitized_task_label = sanitize_label("~{task_identifier}")
        sanitized_task_label = truncate_cluster_name(sanitized_task_label)

        # Generate a unique cluster name
        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'{sanitized_task_label}-hail-{str(uuid.uuid4())[0:13]}'
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
                --worker-machine-type n2-highmem-8
                --master-machine-type n1-highmem-32
                --max-idle=60m --max-age=1440m
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
                        --driver_memory ~{driver_memory} --driver_cores ~{driver_cores}
                        --reference_genome ~{reference_genome}
                        --output_gs_url gs://{cluster_staging_bucket}/{cluster_name}
                        --min_partitions ~{min_partitions}
                        --vcf_url ~{vcf_url}
                        --pca_scores_url ~{pca_scores_url}
                        --min_individual_maf ~{--min_individual_maf}
                        --statistics ~{statistics} --min_kinship ~{min_kinship}
                        --block_size ~{block_size}
                        --task_identifier ~{task_identifier}
                    """

                    print("Running: " + submit_cmd)

                    # Replace newline characters with spaces and remove extra spaces
                    submit_cmd = ' '.join(submit_cmd.split())
                    os.popen(submit_cmd).read()

                    print("Copying results out of staging bucket...")
                    staging_cmd = f'gsutil cp -r gs://{cluster_staging_bucket}/{cluster_name}/~{task_identifier}_relatedness.tsv ./~{task_identifier}_relatedness.tsv'
                    print(staging_cmd)
                    os.popen(staging_cmd).read()

                    staging_cmd = f'gsutil cp -r gs://{cluster_staging_bucket}/{cluster_name}/~{task_identifier}_relatedness_flagged_samples.tsv ./~{task_identifier}_relatedness_flagged_samples.tsv'
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
        File relatedness = "~{task_identifier}_relatedness.tsv"
        File relatedness_flagged_samples = "~{task_identifier}_relatedness_flagged_samples.tsv"
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
