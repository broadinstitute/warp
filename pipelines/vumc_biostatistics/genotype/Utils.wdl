version 1.0

task CopyVcfFile {
  input {
    String input_vcf
    String input_vcf_index

    String? project_id
    String target_bucket
    String genoset
  }

  String new_vcf = "~{target_bucket}/~{genoset}/~{basename(input_vcf)}"
  String new_vcf_index = "~{target_bucket}/~{genoset}/~{basename(input_vcf_index)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} cp ~{input_vcf} \
  ~{input_vcf_index} \
  ~{target_bucket}/~{genoset + "/"}

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_vcf = new_vcf
    String output_vcf_index = new_vcf_index
  }
}

task MoveOrCopyPlinkFile {
  input {
    File source_bed
    File source_bim
    File source_fam

    String action = "mv" #set to cp for copy

    String? project_id
    String target_bucket
  }

  String new_bed = "~{target_bucket}/~{basename(source_bed)}"
  String new_bim = "~{target_bucket}/~{basename(source_bim)}"
  String new_fam = "~{target_bucket}/~{basename(source_fam)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{source_bed} \
  ~{source_bim} \
  ~{source_fam} \
  ~{target_bucket}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    File output_bed = new_bed
    File output_bim = new_bim
    File output_fam = new_fam
  }
}
