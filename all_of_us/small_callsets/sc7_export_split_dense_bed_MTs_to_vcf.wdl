version 1.0

# Workflow: export_split_dense_MTs_to_vcf
# Summary: Given SC4's split dense MatrixTables (one per UCSC BED name), export one VCF per BED name.
# This is the first stage in this pipeline that produces VCF output -- there is no VDS-input equivalent.
#
# Output granularity is per-chromosome by way of the pipeline running per-chromosome end to end (one
# chromosome's phased VCF in at SC1, one chromosome's data throughout). This workflow does not split by
# contig itself.
#
# Indexing is deliberately NOT done here. Per guidance, index the resulting VCFs separately with
# tabix/bcftools as a cheap non-Dataproc task rather than paying for Hail's tabix=True.
#
# STATUS: never run. Written with the failure modes found while bringing SC1/SC2/SC4 up already avoided:
#   - input_split_dense_mts_json_path is String, not File. A File input gets localized by Cromwell to a
#     /mnt/disks/cromwell_root/... path, but the Python gsutil-copies it and needs the gs:// URL. This is
#     exactly what broke SC4's first run.
#   - executor_cores / executor_memory / driver_memory have defaults that fit the machine types pinned
#     below. Values exceeding them (e.g. 8 executor cores on a 4-vCPU worker) are not rejected -- YARN
#     leaves the container request pending forever, which presents as a hung job with no error.
#   - This workflow does not pass --staging_bucket. Every sibling script parses it and never uses it;
#     sc7_export_split_dense_bed_MTs_to_vcf.py correspondingly does not accept it.
#
## Requirements/expectations:
## - **This script CAN ONLY BE RUN in Terra and will NOT work in vanilla cromwell**
## - SC4's split dense MatrixTables plus its json manifest
## - The script is tested and expected to run on the GRCh38 reference genome

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow export_split_dense_MTs_to_vcf {
    input {
        # Analysis Parameters
        String input_split_dense_mts_json_path  # gs:// path to SC4's split dense MTs json manifest.
                                                # String, not File -- see STATUS note above.
        String reference_genome  # Reference genome identifier (e.g., "GRCh38")
        String version_label  # Label to version the run (e.g., "phased_test_v2")
        String task_label  # Label to organize output (e.g., "splitmts2vcf")
        String output_gs_url  # Google Cloud Storage path for the output VCFs
        String output_vcfs_json_file_name  # File name of the output VCFs json manifest, suffix .json
        Boolean overwrite_existing_outputs = false  # If false (default), refuse to overwrite existing outputs

        # Spark sizing -- defaults fit the machine types pinned in the task below
        String executor_cores = "4"      # n1-standard-4 worker has 4 vCPU
        String executor_memory = "10g"   # n1-standard-4 worker has 15 GB; leave headroom for YARN overhead
        String driver_memory = "80g"     # n1-highmem-16 master has 104 GB

        # Cluster Parameters
        # Some of these must be taken from Terra workspace cloud information
        Int num_workers  # Number of workers (per shard) for the Hail cluster
        String gcs_project  # Google project ID for Dataproc
        String gcs_subnetwork_name = 'subnetwork'  # Subnetwork name, set if running in Terra Cromwell
        File submission_script  # sc7_export_split_dense_bed_MTs_to_vcf.py -- the script to run on the cluster
        String region = "us-central1"  # Region for Dataproc, set if running in Terra Cromwell

        # VM Parameters
        String hail_docker = "gcr.io/broad-dsde-methods/aou-auxiliary/hail_dataproc_wdl:0.2.125"  # Docker image with Hail and Google Cloud SDK
    }

    String pipeline_version = "aou_beta"

    call export_split_dense_MTs_to_vcf_task {
        input:
            input_split_dense_mts_json_path = input_split_dense_mts_json_path,
            reference_genome = reference_genome,
            version_label = version_label,
            task_label = task_label,
            output_gs_url = output_gs_url,
            output_vcfs_json_file_name = output_vcfs_json_file_name,
            overwrite_existing_outputs = overwrite_existing_outputs,
            executor_cores = executor_cores,
            executor_memory = executor_memory,
            driver_memory = driver_memory,
            gcs_project = gcs_project,
            num_workers = num_workers,
            gcs_subnetwork_name = gcs_subnetwork_name,
            submission_script = submission_script,
            hail_docker = hail_docker,
            region = region
    }

    output {
        File vcfs_json_file = export_split_dense_MTs_to_vcf_task.vcfs_json_file
    }
}

task export_split_dense_MTs_to_vcf_task {
    input {
        String input_split_dense_mts_json_path
        String reference_genome
        String version_label
        String task_label
        String output_gs_url
        String output_vcfs_json_file_name
        Boolean overwrite_existing_outputs
        String executor_cores
        String executor_memory
        String driver_memory
        File submission_script
        String gcs_project
        String region = "us-central1"
        Int num_workers
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name
        String hail_docker
    }

    RuntimeAttr runtime_default = object {
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

        def sanitize_label(label):
            return re.sub(r'[^a-z0-9-]', '-', label.lower())

        def truncate_cluster_name(name):
            name = name.lower()
            if not name[0].isalpha():
                raise ValueError("Cluster name must start with a letter.")
            truncated = name[:10]
            if truncated[-1] == '-':
                truncated = truncated[:-1] + '0'
            return truncated

        sanitized_version_label = sanitize_label("~{version_label}")
        sanitized_task_label = sanitize_label("~{task_label}")
        sanitized_version_label = truncate_cluster_name(sanitized_version_label)
        sanitized_task_label = truncate_cluster_name(sanitized_task_label)

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'{sanitized_version_label}-{sanitized_task_label}-hail-{str(uuid.uuid4())[0:13]}'
        print(f"Cluster Name: {cluster_name}")
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = f"""
                hailctl dataproc start --num-workers ~{num_workers}
                --region ~{region} --project ~{gcs_project} --service-account {account}
                --worker-machine-type n1-standard-4 --master-machine-type n1-highmem-16 --max-idle=60m --max-age=1440m
                --subnet=projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}
                --enable-component-gateway
                {cluster_name}
            """

            print("Starting cluster...")
            cluster_start_cmd = ' '.join(cluster_start_cmd.split())
            print(os.popen(cluster_start_cmd).read())

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    submit_cmd = f"""
                        gcloud dataproc jobs submit pyspark {script_path}
                        --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account}
                        --driver-log-levels root=WARN --
                        --executor_memory ~{executor_memory} --executor_cores ~{executor_cores}
                        --driver_memory ~{driver_memory} --reference_genome ~{reference_genome}
                        --input_split_dense_mts_json_path ~{input_split_dense_mts_json_path}
                        --output_gs_url ~{output_gs_url}
                        --version_label ~{version_label} --task_label ~{task_label}
                        --output_vcfs_json_file_name ~{output_vcfs_json_file_name}
                        --overwrite_existing_outputs ~{overwrite_existing_outputs}
                    """

                    print("Running: " + submit_cmd)
                    submit_cmd = ' '.join(submit_cmd.split())
                    os.popen(submit_cmd).read()

                    print("Copying results out of staging bucket...")
                    staging_cmd = f'gsutil cp -r ~{output_gs_url}/~{version_label}/~{task_label}/~{version_label}.~{task_label}.~{output_vcfs_json_file_name}  ~{version_label}.~{task_label}.~{output_vcfs_json_file_name}'
                    print(staging_cmd)
                    os.popen(staging_cmd).read()
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

    output {
        File vcfs_json_file = "~{version_label}.~{task_label}.~{output_vcfs_json_file_name}"
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}
