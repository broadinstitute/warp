version 1.0

workflow VUMCHailMatrixExtractSamples {
  input {
    String input_hail_mt_path
    Float input_hail_mt_size #almost equals to the size of corresponding vcf.gz file
    File include_samples

    String reference_genome = "GRCh38"

    String target_prefix

    String? project_id
    String target_gcp_folder
  }

  call HailMatrixExtractSamples {
    input:
      input_hail_mt_path = input_hail_mt_path,
      input_hail_mt_size = input_hail_mt_size,
      reference_genome = reference_genome,
      include_samples = include_samples,
      target_prefix = target_prefix,
      project_id = project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    String output_hail_mt_path = HailMatrixExtractSamples.output_hail_mt_path
  }
}

task HailMatrixExtractSamples {
  input {
    String input_hail_mt_path
    Float input_hail_mt_size
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

  Int disk_size = ceil(input_hail_mt_size / 1024 / 1024 / 1024) + 20
  Int total_memory_gb = memory_gb + 2

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")
  String gcs_output_path = gcs_output_dir + "/" + target_prefix

  command <<<

#https://discuss.hail.is/t/i-get-a-negativearraysizeexception-when-loading-a-plink-file/899
export PYSPARK_SUBMIT_ARGS="--driver-java-options '-XX:hashCode=0' --conf 'spark.executor.extraJavaOptions=-XX:hashCode=0' pyspark-shell"

python3 <<CODE

import hail as hl
import pandas as pd

eligible_samples = pd.read_csv("~{include_samples}", header=None, names=['s'])['s'].tolist()
print(f"There are {len(eligible_samples)} eligible samples.")

hl.init(spark_conf={"spark.driver.memory": "~{memory_gb}g"}, idempotent=True)
hl.default_reference("~{reference_genome}")

mt = hl.read_matrix_table("~{input_hail_mt_path}")

idx = [i for i, s in enumerate(mt.s.collect()) if s in eligible_samples]
mt_eligible = mt.choose_cols(idx)
print(f"Find {mt_eligible.count_cols()} samples in HailMatrix table.")

mt_eligible.write("~{target_prefix}", overwrite=True)

CODE

gsutil ~{"-u " + project_id} -m rsync -Cr ~{target_prefix} ~{gcs_output_path}

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
