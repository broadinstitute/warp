version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow SubsetPhasedVcfsForFlare {
    input {
        # Analysis parameters
        String phased_vcf_gcs_dir          # e.g. gs://prod-drc-broad/aou_phasing/v9
        String chrom                       # e.g. chr21
        String samples_tsv                 # e.g. gs://.../sas_amr.tsv
        String sample_id_column            # e.g. research_id

        String task_identifier             # e.g. sas_amr_chr21_round2
        Int min_partitions = 1200

        # Output locations (GCS paths for the job to write to)
        # NOTE: output_vcf MUST be a gs://... path
        String output_vcf                  # e.g. gs://bucket/path/chr21.subset.vcf.bgz

        # Optional checkpoint output (GCS path, dir ending in .mt)
        Boolean do_checkpoint = false
        String? checkpoint_path            # e.g. gs://bucket/path/mt_sub_chr21.ckpt.mt
        Boolean checkpoint_overwrite = true

        # Entry fields
        String entry_fields = "GT"         # comma-separated; "" means keep all
        Boolean force_bgz = false          # passed through to hl.import_vcf

        # Hail/Spark resources (passed to your script)
        String executor_cores
        String driver_cores
        String executor_memory
        String driver_memory
        String reference_genome = "GRCh38"

        # Cluster params
        Int num_workers
        String gcs_project
        String gcs_subnetwork_name = "subnetwork"
        String region = "us-central1"
        Int max_idle = 60
        Int max_age = 1440

        # Script & docker
        File submission_script
        String hail_docker = "us.gcr.io/broad-dsde-methods/lichtens/hail_dataproc_wdl:1.1"

        RuntimeAttr? runtime_attr_override
    }
	String pipeline_version = "aou_10.0.0"

    call subset_phased_vcf_task {
        input:
            phased_vcf_gcs_dir = phased_vcf_gcs_dir,
            chrom = chrom,
            samples_tsv = samples_tsv,
            sample_id_column = sample_id_column,
            task_identifier = task_identifier,
            min_partitions = min_partitions,
            output_vcf = output_vcf,
            do_checkpoint = do_checkpoint,
            checkpoint_path = checkpoint_path,
            checkpoint_overwrite = checkpoint_overwrite,
            entry_fields = entry_fields,
            force_bgz = force_bgz,
            executor_cores = executor_cores,
            driver_cores = driver_cores,
            executor_memory = executor_memory,
            driver_memory = driver_memory,
            reference_genome = reference_genome,
            gcs_project = gcs_project,
            num_workers = num_workers,
            gcs_subnetwork_name = gcs_subnetwork_name,
            submission_script = submission_script,
            hail_docker = hail_docker,
            region = region,
            max_idle = max_idle,
            max_age = max_age,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File subset_vcf = subset_phased_vcf_task.subset_vcf
        File subset_vcf_tbi = subset_phased_vcf_task.subset_vcf_tbi
    }
}

task subset_phased_vcf_task {
    input {
        String phased_vcf_gcs_dir
        String chrom
        String samples_tsv
        String sample_id_column
        String task_identifier
        Int min_partitions
        String output_vcf

        Boolean do_checkpoint
        String? checkpoint_path
        Boolean checkpoint_overwrite

        String entry_fields
        Boolean force_bgz

        String executor_cores
        String executor_memory
        String driver_cores
        String driver_memory
        String reference_genome

        File submission_script
        String gcs_project
        String region
        Int num_workers
        String gcs_subnetwork_name
        String hail_docker
        Int max_idle
        Int max_age

        RuntimeAttr? runtime_attr_override
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

		# Comfigure gcloud account
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
		import sys
		import uuid
		import re
		from google.cloud import dataproc_v1 as dataproc
		from datetime import datetime

		def popen_read_checked(cmd):
			with os.popen(cmd) as stream:
				output = stream.read()
				status = stream.close()  # returns None on success, exit code << 8 on failure
				if status is not None:   # means command failed
					raise RuntimeError(f"Command failed with exit code {status >> 8}: {cmd}\n{output}")
				return output

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
				--worker-machine-type n1-standard-4
				--master-machine-type n1-highmem-32
				--max-idle=~{max_idle}m --max-age=~{max_age}m
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
					# Submit the job
					submit_cmd = f"""
						gcloud dataproc jobs submit pyspark {script_path}
						--cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account}
						--driver-log-levels root=WARN --
						--executor_memory ~{executor_memory} --executor_cores ~{executor_cores}
						--driver_memory ~{driver_memory} --driver_cores ~{driver_cores}
						--reference_genome ~{reference_genome}
						--min_partitions ~{min_partitions}
						--tmp_dir {tmp_dir}
						--phased_vcf_gcs_dir ~{phased_vcf_gcs_dir}
						--chrom ~{chrom}
						--samples_tsv ~{samples_tsv}
						--sample_id_column ~{sample_id_column}
						--entry_fields "~{entry_fields}"
						{"--force_bgz" if "~{force_bgz}" == "true" else ""}
						{ckpt_args}
						--output_vcf ~{output_vcf}
						--log_level INFO
					"""
					print("Running: " + submit_cmd)
					submit_cmd = " ".join(submit_cmd.split())
					print("Running:", submit_cmd)
					print(popen_read_checked(submit_cmd))

					# Copy outputs back locally for WDL outputs
					# output_vcf is a gs:// path; Cromwell needs local files as task outputs
					out_vcf = "~{task_identifier}.subset.vcf.bgz"
					out_tbi = out_vcf + ".tbi"

					copy_vcf_cmd = f"gsutil -m cp ~{output_vcf} ./{out_vcf}"
					copy_tbi_cmd = f"gsutil -m cp ~{output_vcf}.tbi ./{out_tbi}"
					print(popen_read_checked(copy_vcf_cmd))
					print(popen_read_checked(copy_tbi_cmd))

					break

		except Exception as e:
			timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
			print(timestamp, "Exception raised!")
			print(e)
			raise
		finally:
			print(f"Stopping cluster: {cluster_name}")
			os.popen(
				"gcloud dataproc clusters delete --quiet --project {} --region {} --account {} {}".format(
					"~{gcs_project}", "~{region}", account, cluster_name
				)
			).read()
		EOF

		echo "Complete"
	>>>

    output {
        File subset_vcf = "~{task_identifier}.subset.vcf.bgz"
        File subset_vcf_tbi = "~{task_identifier}.subset.vcf.bgz.tbi"
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
