version 1.0

workflow VUMCHailMTIncludeSamples {
  input {
    String input_hail_mt_path
    File include_samples

    String reference_genome = "GRCh38"

    String target_prefix

    String? project_id
    String target_gcp_folder
  }

  call HailMTIncludeSamples {
    input:
      input_hail_mt_path = input_hail_mt_path,
      reference_genome = reference_genome,
      include_samples = include_samples,
      target_prefix = target_prefix,
      project_id = project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    String output_hail_mt_path = HailMTIncludeSamples.output_hail_mt_path
  }
}

task HailMTIncludeSamples {
  input {
    String input_hail_mt_path
    String reference_genome
    File include_samples
    String target_prefix
    String? project_id
    String target_gcp_folder

    String docker = "shengqh/hail_gcp:20240211"
    Int memory_gb = 64
    Int preemptible = 1
    Int cpu = 4
    Int boot_disk_gb = 10  
  }

  Int disk_size = 20
  Int total_memory_gb = memory_gb + 2

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")
  String gcs_output_path = gcs_output_dir + "/" + target_prefix

  command <<<

#https://discuss.hail.is/t/i-get-a-negativearraysizeexception-when-loading-a-plink-file/899
export PYSPARK_SUBMIT_ARGS="--driver-java-options '-XX:hashCode=0' --conf 'spark.executor.extraJavaOptions=-XX:hashCode=0' pyspark-shell"

python3 <<CODE

import hail as hl

hl.init(spark_conf={"spark.driver.memory": "~{memory_gb}g"}, default_reference="~{reference_genome}", idempotent=True)

mt = hl.read_matrix_table("~{input_hail_mt_path}")
single_mt = mt.filter_cols(mt.s == 'R200013600')
single_mt.write("~{gcs_output_path}", overwrite=True)

CODE

>>>

  runtime {
    cpu: cpu
    docker: "~{docker}"
    preemptible: preemptible
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{total_memory_gb} GiB"
    bootDiskSizeGb: boot_disk_gb
  }
  output {
    String output_hail_mt_path = "~{gcs_output_path}"
  }
}
