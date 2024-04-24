version 1.0

workflow VUMCVcf2HailMatrix {
  #modified based on https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/ConvertToHailMTT2T:hangsu_phasing?tab=files
 
  input {
    File source_vcf
    File source_vcf_index

    String reference_genome = "GRCh38"

    String target_prefix

    String? project_id
    String target_gcp_folder
  }

  call Vcf2HailMatrix {
    input:
      source_vcf = source_vcf,
      source_vcf_index = source_vcf_index,
      reference_genome = reference_genome,
      target_prefix = target_prefix,
      project_id = project_id,
      target_gcp_folder = target_gcp_folder,
  }

  output {
    String hail_gcs_path = Vcf2HailMatrix.hail_gcs_path
  }
}

task Vcf2HailMatrix {
  meta {
    description: "Convert a .vcf.bgz file to a Hail MatrixTable and copy it to a final gs:// URL."
  }

  parameter_meta {
    source_vcf: "The input .vcf.bgz file."
    source_vcf_index: "The input .vcf.bgz.tbi file."
    reference_genome: "The reference genome to use.  Currently only GRCh38 is supported."
    target_prefix: "The prefix to use for the output MatrixTable."

    project_id: "The GCP project to use for gsutil."
    target_gcp_folder: "The output GCP directory to copy the MatrixTable to."

    docker: "The docker image to use for this task."
    memory_gb: "The amount of memory to use for this task."
    preemptible: "Number of preemptible tries to use for this task."
  }  

  input {
    File source_vcf
    File source_vcf_index
    String reference_genome
    String target_prefix

    String? project_id
    String target_gcp_folder

    String docker = "shengqh/hail_gcp:20240211"
    Int memory_gb = 64
    Int preemptible = 1
    Int cpu = 4
    Int boot_disk_gb = 10
  }

  Int disk_size = 100 + 3*ceil(size(source_vcf, "GB"))
  Int total_memory_gb = memory_gb + 2

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")
  String gcs_output_path = gcs_output_dir + "/" + target_prefix

  command <<<

#https://discuss.hail.is/t/i-get-a-negativearraysizeexception-when-loading-a-plink-file/899
export PYSPARK_SUBMIT_ARGS="--driver-java-options '-XX:hashCode=0' --conf 'spark.executor.extraJavaOptions=-XX:hashCode=0' pyspark-shell"

python3 <<CODE

import hail as hl

hl.init(spark_conf={"spark.driver.memory": "~{memory_gb}g"})

callset = hl.import_vcf("~{source_vcf}",
                        array_elements_required=False,
                        force_bgz=True,
                        reference_genome='~{reference_genome}')

callset.write("~{target_prefix}", overwrite=True)

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
    String hail_gcs_path = "~{gcs_output_path}"
  }
}
