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
