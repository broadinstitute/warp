version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils

workflow VUMCRegenie {
  input {
    File pgen_file
    File pvar_file
    File psam_file

    File phenoFile
    String phenoColList

    File covarFile
    String covarColList

    String target_prefix

    String qc_option="--mac 100 --geno 0.1 --maf 0.1 --max-maf 0.9 --hwe 1e-15 --snps-only --not-chr 23-27"

    String? project_id
    String? target_gcp_folder
  }

  call PgenQCFilter {
    input:
      input_pgen = pgen_file,
      input_pvar = pvar_file,
      input_psam = psam_file,
      target_prefix = target_prefix,
      qc_option = qc_option
  }

  call Regenie {
    input:
      input_qc_pgen = PgenQCFilter.output_pgen,
      input_qc_pvar = PgenQCFilter.output_pvar,
      input_qc_psam = PgenQCFilter.output_psam,
      input_pgen = pgen_file,
      input_pvar = pvar_file,
      input_psam = psam_file,
      phenoFile = phenoFile,
      phenoColList = phenoColList,
      covarFile = covarFile,
      covarColList = covarColList,
      target_prefix = target_prefix
  }

  output {
    File pred_list = Regenie.pred_list
    Array[File] loco_files = Regenie.loco_files
    Array[File] regenie_files = Regenie.regenie_files
  }
}

task PgenQCFilter {
  input {
    File input_pgen
    File input_pvar
    File input_psam

    String target_prefix

    String qc_option

    Int memory_gb = 20
    Int cpu = 8

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([input_pgen, input_pvar, input_psam], "GB")  * 2) + 20

  String new_pgen = target_prefix + ".qc.pgen"
  String new_pvar = target_prefix + ".qc.pvar"
  String new_psam = target_prefix + ".qc.psam"

  command <<<

plink2 \
  --pgen ~{input_pgen} \
  --pvar ~{input_pvar} \
  --psam ~{input_psam} \
  ~{qc_option} \
  --make-pgen \
  --out ~{target_prefix}.qc

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen = new_pgen
    File output_pvar = new_pvar
    File output_psam = new_psam
  }
}

task Regenie {
  input {
    File input_qc_pgen
    File input_qc_pvar
    File input_qc_psam

    File input_pgen
    File input_pvar
    File input_psam

    File phenoFile
    String phenoColList

    File covarFile
    String covarColList

    String step1_option = "--loocv --bsize 1000 --lowmem"
    String step2_option = "--minINFO 0.3 --minMAC 30 --bsize 1000"

    String target_prefix

    Int memory_gb = 20
    Int cpu = 8

    String docker = "skoyamamd/regenie:3.4.2"
  }

  Int disk_size = ceil(size([input_qc_pgen, input_qc_pvar, input_qc_psam, input_pgen, input_pvar, input_psam], "GB")) + 20

  String qc_pgen_prefix = sub(input_qc_pgen, ".pgen$", "")
  String pgen_prefix = sub(input_pgen, ".pgen$", "")

  command <<<

regenie --step 1 \
  --qt \
  --pgen ~{qc_pgen_prefix} \
  -p ~{phenoFile} \
  --phenoColList ~{phenoColList} \
  -c ~{covarFile} \
  --covarColList ~{covarColList} \
  ~{step1_option} \
  --threads ~{cpu} \
  --out ~{target_prefix}.step1 \
  --force-step1
  
regenie --step 2 \
  --qt \
  --pgen ~{pgen_prefix} \
  -p ~{phenoFile} \
  --phenoColList ~{phenoColList} \
  -c ~{covarFile} \
  --covarColList ~{covarColList} \
  ~{step2_option} \
  --threads ~{cpu} \
  --pred "~{target_prefix}.step1_pred.list" \
  --out ~{target_prefix}.step2

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File pred_list = "~{target_prefix}.step1_pred.list"
    Array[File] loco_files = glob("~{target_prefix}.step1*.loco")
    Array[File] regenie_files = glob("~{target_prefix}.step2*.regenie")
  }
}
