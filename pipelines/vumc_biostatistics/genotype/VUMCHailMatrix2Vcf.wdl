version 1.0

workflow VUMCHailMatrix2Vcf {
  input {
    String input_hail_mt_path
    Float expect_vcf_size #almost equals to the size of corresponding vcf.gz file

    String reference_genome = "GRCh38"

    String target_prefix

    String? project_id
    String target_gcp_folder
  }

  call HailMatrix2Vcf {
    input:
      input_hail_mt_path = input_hail_mt_path,
      expect_vcf_size = expect_vcf_size,
      reference_genome = reference_genome,
      target_prefix = target_prefix,
      project_id = project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    String output_vcf = HailMatrix2Vcf.output_vcf
    String output_vcf_tbi = HailMatrix2Vcf.output_vcf_tbi
  }
}

task HailMatrix2Vcf {
  input {
    String input_hail_mt_path
    Float expect_vcf_size
    String reference_genome
    String target_prefix
    String? project_id
    String target_gcp_folder

    String docker = "shengqh/hail_gcp:20240213"
    Int memory_gb = 64
    Int preemptible = 1
    Int cpu = 4
    Int boot_disk_gb = 10  
    Float disk_size_factor = 3
  }

  Int disk_size = ceil(expect_vcf_size / 1024 / 1024 / 1024 * disk_size_factor) + 20
  Int total_memory_gb = memory_gb + 2

  #for hail export, we have to use .bgz as the output suffix
  String hail_vcf = target_prefix + ".vcf.bgz"

  #in order to be consistent with the original VCF file, we use .vcf.gz as the output suffix
  String target_vcf = target_prefix + ".vcf.gz"

  String target_vcf_tbi = target_vcf + ".tbi"
  String target_sample_file = target_vcf + ".samples.txt"
  
  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")
  String gcs_output_vcf = gcs_output_dir + "/" + target_vcf
  String gcs_output_vcf_tbi = gcs_output_vcf + ".tbi"
  String gcs_output_sample_file = gcs_output_vcf + ".samples.txt"

  command <<<

#https://discuss.hail.is/t/i-get-a-negativearraysizeexception-when-loading-a-plink-file/899
export PYSPARK_SUBMIT_ARGS="--driver-java-options '-XX:hashCode=0' --conf 'spark.executor.extraJavaOptions=-XX:hashCode=0' pyspark-shell"

mkdir -p ./tmp

python3 <<CODE

import hail as hl

hl.init(spark_conf={
    "spark.driver.memory": "~{memory_gb}g",
    "spark.local.dir": "./tmp"
  },
  tmp_dir="./tmp",
  local_tmpdir="./tmp",
  idempotent=True)
hl.default_reference("~{reference_genome}")

mt = hl.read_matrix_table("~{input_hail_mt_path}")
hl.export_vcf(mt, "~{hail_vcf}", tabix = False)

CODE

mv ~{hail_vcf} ~{target_vcf}

#bcftools cannot query number of records based on hail index file.
#We use bcftools index to generate tbi file
bcftools index -t ~{target_vcf} --threads ~{cpu} -o ~{target_vcf_tbi}

bcftools query -l ~{target_vcf} > ~{target_sample_file}

cat ~{target_sample_file} | wc -l > num_samples.txt

bcftools index -n ~{target_vcf} > num_variants.txt

gsutil ~{"-u " + project_id} -m cp ~{target_vcf} ~{gcs_output_vcf}

gsutil ~{"-u " + project_id} -m cp ~{target_vcf_tbi} ~{gcs_output_vcf_tbi}

gsutil ~{"-u " + project_id} -m cp ~{target_sample_file} ~{gcs_output_sample_file}

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
    String output_vcf = "~{gcs_output_vcf}"
    String output_vcf_tbi = "~{gcs_output_vcf_tbi}"
    String output_vcf_sample = "~{gcs_output_sample_file}"
    Int output_vcf_num_samples = read_int("num_samples.txt")
    Int output_vcf_num_variants = read_int("num_variants.txt")
  }
}
