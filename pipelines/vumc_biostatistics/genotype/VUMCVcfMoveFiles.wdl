version 1.0

workflow VUMCVcfMoveFiles {
  input {
    String input_vcf
    String input_vcf_index
    String input_vcf_sample

    String? project_id
    String target_gcp_folder
  }

  call VcfMoveFiles {
    input:
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      input_vcf_sample = input_vcf_sample,
      project_id = project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    File output_vcf = VcfMoveFiles.output_vcf
    File output_vcf_index = VcfMoveFiles.output_vcf_index
    File output_vcf_sample = VcfMoveFiles.output_vcf_sample
  }
}

task VcfMoveFiles {
  input {
    String input_vcf
    String input_vcf_index
    String input_vcf_sample

    String? project_id
    String target_gcp_folder
  }

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")

  String gcs_output_vcf = "~{gcs_output_dir}/~{basename(input_vcf)}"
  String gcs_output_vcf_index = "~{gcs_output_dir}/~{basename(input_vcf_index)}"
  String gcs_output_sample_file ="~{gcs_output_dir}/~{basename(input_vcf_sample)}"

  command <<<

echo gsutil ~{"-u " + project_id} -m mv ~{input_vcf} ~{input_vcf_index} ~{input_vcf_sample} "~{gcs_output_dir}/"
gsutil ~{"-u " + project_id} -m mv ~{input_vcf} ~{input_vcf_index} ~{input_vcf_sample} "~{gcs_output_dir}/"

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    # The output has to be defined as String, otherwise the file would be delocalized and failed.
    String output_vcf = "~{gcs_output_vcf}"
    String output_vcf_index = "~{gcs_output_vcf_index}"
    String output_vcf_sample = "~{gcs_output_sample_file}"
  }
}
