version 1.0

task MoveOrCopyVcfFile {
  input {
    String input_vcf
    String input_vcf_index

    Boolean is_move_file = true

    String? project_id
    String target_bucket
    String genoset
    String? GRID
  }

  String action = if (is_move_file) then "mv" else "cp"

  String target_folder = if(defined(GRID)) then "~{target_bucket}/~{genoset}/~{GRID}" else "~{target_bucket}/~{genoset}"
  String new_vcf = "~{target_folder}/~{basename(input_vcf)}"
  String new_vcf_index = "~{target_folder}/~{basename(input_vcf_index)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{input_vcf} \
  ~{input_vcf_index} \
  ~{target_folder}/

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
    String source_bed
    String source_bim
    String source_fam

    Boolean is_move_file = true

    String? project_id
    String target_bucket
  }

  String action = if (is_move_file) then "mv" else "cp"

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
    String output_bed = new_bed
    String output_bim = new_bim
    String output_fam = new_fam
  }
}
