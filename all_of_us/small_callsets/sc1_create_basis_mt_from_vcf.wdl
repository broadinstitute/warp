version 1.0

# Workflow: filter_vcf_to_MT_by_bed
# Summary: Given one or more phased, per-chromosome VCFs and a list of UCSC BED files, this workflow
# creates a MatrixTable (MT) of filtered variants corresponding to regions in these UCSC BED files.
#
# This is the phased-VCF-input counterpart to sc1_create_basis_mt_from_vds.wdl -- mirrors its structure
# (Dataproc cluster launch pattern, runtime defaults) exactly, swapping vds_url for vcf_url and adding the
# inputs specific to sc1_create_basis_mt_from_vcf.py (test_intervals, sample_ids,
# overwrite_existing_outputs). See that script/WDL for the VDS-input sibling.
#
# STUB STATUS: structurally complete and mirrors the validated sc1-VDS WDL, but not yet run through
# Cromwell/Terra -- needs an actual end-to-end WDL execution test before being considered validated,
# separate from the Python script's own notebook-based validation (see
# sc1_create_basis_mt_from_vcf_validation_notes.md).
#
## Requirements/expectations:
## - **This script CAN ONLY BE RUN in Terra and will NOT work in vanilla cromwell**
## - Phased, per-chromosome VCFs (bgzipped), one chromosome per run
## - UCSC BED files to specify genomic regions of interest: Exome, ClinVar, ACAF Threshold
## - The script is tested and expected to run on the GRCh38 reference genome

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filter_vcf_to_MT_by_bed {
    input {
        # Analysis Parameters
        String vcf_url  # Path/glob, or comma/space-delimited list of paths, to the input phased VCF(s)
        File exome_bed_path  # Path to the exome UCSC BED file
        File clinvar_bed_path  # Path to the ClinVar UCSC BED file
        File acaf_threshold_bed_path  # Path to the aCAF threshold UCSC BED file
        String reference_genome  # Reference genome identifier (e.g., "GRCh38")
        Boolean is_test  # Flag indicating if this is a test run
        String test_chrs = ""  # Chromosome(s) to filter to, only used if vcf_url spans multiple chromosomes
        String test_intervals = ""  # Comma-delimited locus interval(s) for a fast smoke test
        String sample_ids = ""  # Comma/space-delimited sample IDs to subset to, for a fast test
        String version_label  # Label to version the run (e.g., "aou_phased_delta")
        String task_label  # Label to organize output (e.g., "vcf2mt")
        Int min_partitions  # Minimum number of partitions for optimization
        Boolean is_filtering_FT  # Whether to filter GTs based on "Failed" FT field, if present
        Boolean is_keep_as_vqsr_fields  # Flag to keep VQSR fields, if present
        Int? max_alt_alleles  # Maximum number of alternate alleles for filtering
        String executor_cores  # Number of cores assigned to each Spark executor
        String executor_memory  # Memory assigned to each Spark executor
        String driver_memory  # Memory assigned to the Spark driver
        String output_gs_url  # Google Cloud Storage path for the output MatrixTable
        String basis_mt_json_file_name  # File name of the basis dense MatrixTables json file, suffix should be .json
        String fields_to_drop  # Fields to drop from the VCF-derived MT
        Boolean overwrite_existing_outputs = false  # If false (default), refuse to overwrite existing outputs

        # Cluster Parameters
        # Some of these must be taken from Terra workspace cloud information
        Int num_workers  # Number of workers (per shard) for the Hail cluster
        String gcs_project  # Google project ID for Dataproc
        String gcs_subnetwork_name = 'subnetwork'  # Subnetwork name, set if running in Terra Cromwell
        File submission_script  # sc1_create_basis_mt_from_vcf.py -- the script to run on the cluster
        String region = "us-central1"  # Region for Dataproc, set if running in Terra Cromwell

        # VM Parameters
        String hail_docker = "gcr.io/broad-dsde-methods/aou-auxiliary/hail_dataproc_wdl:0.2.125"  # Docker image with Hail and Google Cloud SDK
    }

    String pipeline_version = "aou_beta"

    call filter_vcf_and_export_as_MT {
        # Task inputs mirror workflow inputs
        input:
            vcf_url = vcf_url,
            exome_bed_path = exome_bed_path,
            clinvar_bed_path = clinvar_bed_path,
            acaf_threshold_bed_path = acaf_threshold_bed_path,
            reference_genome = reference_genome,
            is_test = is_test,
            test_chrs = test_chrs,
            test_intervals = test_intervals,
            sample_ids = sample_ids,
            version_label = version_label,
            task_label = task_label,
            min_partitions = min_partitions,
            is_filtering_FT = is_filtering_FT,
            is_keep_as_vqsr_fields = is_keep_as_vqsr_fields,
            max_alt_alleles = max_alt_alleles,
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
            basis_mt_json_file_name = basis_mt_json_file_name,
            fields_to_drop = fields_to_drop,
            overwrite_existing_outputs = overwrite_existing_outputs
    }

    output {
        File basis_mt_json_file = filter_vcf_and_export_as_MT.basis_mt_json_file
    }
}

task filter_vcf_and_export_as_MT {
    input {
        # Task-specific inputs with descriptions
        String vcf_url
        File exome_bed_path
        File clinvar_bed_path
        File acaf_threshold_bed_path
        String reference_genome
        Boolean is_test
        String test_chrs
        String test_intervals
        String sample_ids
        String version_label
        String task_label
        Int min_partitions
        Boolean is_filtering_FT
        Boolean is_keep_as_vqsr_fields
        Int? max_alt_alleles
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
        String basis_mt_json_file_name
        String fields_to_drop
        Boolean overwrite_existing_outputs
    }

    parameter_meta {
        # Additional metadata for input parameters
        # no need to localize gs url files as Hail can read from GCS
        exome_bed_path: {localization_optional: true}
        clinvar_bed_path: {localization_optional: true}
        acaf_threshold_bed_path: {localization_optional: true}
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
                        --exome_bed_path ~{exome_bed_path} --clinvar_bed_path ~{clinvar_bed_path}
                        --acaf_threshold_bed_path ~{acaf_threshold_bed_path} --is_test ~{is_test}
                        --output_gs_url ~{output_gs_url}
                        --basis_mt_json_file_name ~{basis_mt_json_file_name}
                        --test_chrs '~{test_chrs}' --test_intervals '~{test_intervals}'
                        --sample_ids '~{sample_ids}'
                        --version_label ~{version_label}
                        --min_partitions ~{min_partitions} --task_label ~{task_label}
                        --vcf_url ~{vcf_url} --is_filtering_FT ~{is_filtering_FT} --is_keep_as_vqsr_fields ~{is_keep_as_vqsr_fields}
                        ~{"--max_alt_alleles " + max_alt_alleles}
                        --fields_to_drop '~{fields_to_drop}'
                        --overwrite_existing_outputs ~{overwrite_existing_outputs}
                        --staging_bucket gs://{cluster_staging_bucket}/{cluster_name}
                    """

                    print("Running: " + submit_cmd)

                    # Replace newline characters with spaces and remove extra spaces
                    submit_cmd = ' '.join(submit_cmd.split())
                    os.popen(submit_cmd).read()

                    print("Copying results out of staging bucket...")
                    staging_cmd = f'gsutil cp -r ~{output_gs_url}/~{version_label}.~{task_label}.~{basis_mt_json_file_name}  ~{version_label}.~{task_label}.~{basis_mt_json_file_name}'
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
        File basis_mt_json_file = "~{version_label}.~{task_label}.~{basis_mt_json_file_name}"
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
