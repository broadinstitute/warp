version 1.0

import "CreateVcfIndex.wdl"

struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
}

workflow RunRemovePhasedSamples {
  meta {
    description: "Scatter Dataproc jobs over input MT paths and run remove_phased_samples.py on each shard, with optional chrX mode."
    allowNestedInputs: true
  }

  input {
    # Analysis inputs
    Array[String] input_mt_paths
    String remove_samples_tsv
    String remove_id_col = "research_id"
    Boolean run_chrX = false
    String? participant_sex_tsv
    String? metadata_vcf_or_header
    Boolean write_out_mt = false
    Boolean overwrite = false

    # Cluster orchestration
    String gcs_project
    String gcs_subnetwork_name = "subnetwork"
    String region = "us-central1"
    String output_bucket_path
    File submission_script

    # Spark executor/driver config (passed as script args to remove_phased_samples.py)
    Int executor_cores = 4
    String executor_memory = "26g"
    Int driver_cores = 32
    String driver_memory = "60g"
    Int spark_task_max_failures = 20

    # WDL runtime for the lightweight launcher VM
    String hail_docker = "gcr.io/broad-dsde-methods/aou-auxiliary/hail_dataproc_wdl:0.2.134"
  }

  String pipeline_version = "aou_9.0.2"
  String output_bucket_path_with_trailing_slash = sub(output_bucket_path, "/$", "") + "/"

  scatter (path in input_mt_paths) {
    call RemovePhasedSamplesOnDataproc {
      input:
        input_mt_path = path,
        remove_samples_tsv = remove_samples_tsv,
        remove_id_col = remove_id_col,
        run_chrX = run_chrX,
        participant_sex_tsv = participant_sex_tsv,
        metadata_vcf_or_header = metadata_vcf_or_header,
        write_out_mt = write_out_mt,
        overwrite = overwrite,
        gcs_project = gcs_project,
        gcs_subnetwork_name = gcs_subnetwork_name,
        region = region,
        output_bucket = output_bucket_path_with_trailing_slash,
        submission_script = submission_script,
        executor_cores = executor_cores,
        executor_memory = executor_memory,
        driver_cores = driver_cores,
        driver_memory = driver_memory,
        spark_task_max_failures = spark_task_max_failures,
        hail_docker = hail_docker
    }

    call CreateVcfIndex.CreateVcfIndex as IndexVcf {
      input:
        vcf_input = RemovePhasedSamplesOnDataproc.filtered_vcf_url
    }
  }

  output {
    Array[String] filtered_vcf_urls = RemovePhasedSamplesOnDataproc.filtered_vcf_url
    Array[String] filtered_vcf_index_urls = IndexVcf.output_vcf_index
    Array[String?] filtered_mt_urls = RemovePhasedSamplesOnDataproc.filtered_mt_url
  }
}

task RemovePhasedSamplesOnDataproc {
  input {
    String input_mt_path
    String remove_samples_tsv
    String remove_id_col
    Boolean run_chrX = false
    String? participant_sex_tsv
    String? metadata_vcf_or_header
    Boolean write_out_mt
    Boolean overwrite

    File submission_script
    String output_bucket

    String gcs_project
    String gcs_subnetwork_name
    String region = "us-central1"

    String master_machine_type = "n1-highmem-32"
    Float master_memory_fraction = 0.8
    String worker_machine_type = "n1-highmem-8"
    Int num_workers = 16
    Int num_preemptible_workers = 0
    Int time_to_live_minutes = 14400
    RuntimeAttr? runtime_attr_override

    # Spark executor/driver config (passed as script args to remove_phased_samples.py)
    Int executor_cores = 4
    String executor_memory = "26g"
    Int driver_cores = 32
    String driver_memory = "60g"
    Int spark_task_max_failures = 8

    String hail_docker
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 6.5,
    disk_gb: 15,
    cpu_cores: 1,
    preemptible_tries: 0,
    max_retries: 0,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  String mt_path_no_trailing_slash = sub(input_mt_path, "/$", "")
  String mt_basename = sub(mt_path_no_trailing_slash, "^.*/", "")
  String mt_output_stem = sub(mt_basename, "\\.mt$", "")

  String output_filtered_vcf_url = output_bucket + mt_output_stem + ".filtered.vcf.bgz"
  String output_filtered_mt_url = output_bucket + mt_output_stem + ".filtered.mt"

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
    print("Running remove_phased_samples python code...")
    import os
    import uuid
    from google.cloud import dataproc_v1 as dataproc

    run_chrx = "~{run_chrX}".lower() == "true"
    participant_sex_tsv = "~{if defined(participant_sex_tsv) then select_first([participant_sex_tsv]) else ""}"

    if run_chrx and not participant_sex_tsv:
        raise Exception("run_chrX=true requires participant_sex_tsv to be provided.")

    sex_tsv_arg = f"--sex-tsv {participant_sex_tsv}" if run_chrx else ""

    cluster_prefix = "~{mt_output_stem}".lower().replace("_", "-").replace(".", "-")
    cluster_prefix = cluster_prefix[:20] if len(cluster_prefix) > 20 else cluster_prefix
    cluster_name = f"rmph-{cluster_prefix}-hail-step1-{str(uuid.uuid4())[0:13]}"

    script_path = "~{submission_script}"

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
            raise Exception("Failed to start cluster successfully")

        cluster_client = dataproc.ClusterControllerClient(
            client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
        )

        for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
            if cluster.cluster_name == cluster_name:
                cluster_temp_bucket = cluster.config.temp_bucket

                submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} \
                --driver-log-levels root=WARN \
                -- \
                --mt-path ~{input_mt_path} \
                --remove-samples-tsv ~{remove_samples_tsv} \
                --remove-id-col ~{remove_id_col} \
                {sex_tsv_arg} \
                --out-vcf ~{output_filtered_vcf_url} \
                --spark-executor-cores ~{executor_cores} \
                --spark-executor-memory ~{executor_memory} \
                --spark-driver-cores ~{driver_cores} \
                --spark-driver-memory ~{driver_memory} \
                --spark-task-max-failures ~{spark_task_max_failures} \
                ~{if write_out_mt then "--out-mt-path " + output_filtered_mt_url else ""} \
                ~{if defined(metadata_vcf_or_header) then "--metadata-vcf-or-header " + select_first([metadata_vcf_or_header]) else ""} \
                ~{if overwrite then "--overwrite" else ""} \
                --temp-bucket gs://{cluster_temp_bucket}/{cluster_name}'''

                print("Running: " + submit_cmd)
                f = os.popen(submit_cmd)
                f.read()
                if (f.close() != None):
                    raise Exception("Failed to submit cluster job successfully")

                break

    except Exception as e:
        print(e)
        raise
    finally:
        print(f"Stopping cluster: {cluster_name}")
        os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} -q {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

    EOF

    echo "Complete"
  >>>

  output {
    String filtered_vcf_url = output_filtered_vcf_url
    String? filtered_mt_url = output_filtered_mt_url
  }

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: hail_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
}
