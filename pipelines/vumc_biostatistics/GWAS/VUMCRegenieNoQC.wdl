version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils

workflow VUMCRegenieNoQC {
  input {
    File qc_filtered_pgen
    File qc_filtered_pvar
    File qc_filtered_psam

    Array[String] chromosomes
    Array[File] pgen_files
    Array[File] pvar_files
    Array[File] psam_files

    File phenoFile
    String phenoColList
    Boolean is_binary_traits = false

    File covarFile
    String covarColList

    String output_prefix

    String? project_id
    String? target_gcp_folder
  }

  call RegenieStep1FitModel {
    input:
      input_qc_pgen = qc_filtered_pgen,
      input_qc_pvar = qc_filtered_pvar,
      input_qc_psam = qc_filtered_psam,
      phenoFile = phenoFile,
      phenoColList = phenoColList,
      is_binary_traits = is_binary_traits,
      covarFile = covarFile,
      covarColList = covarColList,
      output_prefix = output_prefix
  }

  output {
    File model_list_file = RegenieStep1FitModel.model_list_file
    Array[File] loco_files = RegenieStep1FitModel.loco_files
  }
}

task RegenieStep1FitModel {
  input {
    File input_qc_pgen
    File input_qc_pvar
    File input_qc_psam

    File phenoFile
    String phenoColList
    Boolean is_binary_traits

    File covarFile
    String covarColList

    String step1_option = "--bsize 1000 --lowmem"

    String output_prefix

    Int memory_gb = 100
    Int cpu = 8

    String docker = "skoyamamd/regenie:3.4.2"
  }

  Int disk_size = ceil(size([input_qc_pgen, input_qc_pvar, input_qc_psam], "GB")) + 20

  String call_type = if(is_binary_traits) then "--bt" else "--qt"

  command <<<

qc_pgen='~{input_qc_pgen}'
qc_pgen_prefix=${qc_pgen%.*}

regenie --step 1 \
  ~{call_type} \
  --pgen ${qc_pgen_prefix} \
  -p ~{phenoFile} \
  --phenoColList ~{phenoColList} \
  -c ~{covarFile} \
  --covarColList ~{covarColList} \
  ~{step1_option} \
  --threads ~{cpu} \
  --out ~{output_prefix} \
  --force-step1

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File model_list_file = "~{output_prefix}_pred.list"
    Array[File] loco_files=glob("*.loco.gz")    
  }
}
