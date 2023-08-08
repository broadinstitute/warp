version 1.0

workflow VUMCMoveHaplotypecallerResult {
  input {
    String genoset
    String GRID
    String? project_id

    String output_vcf
    String output_vcf_index

    String target_bucket
  }

  call MoveHaplotypecallerResult as mh {
    input:
      genoset = genoset,
      GRID = GRID,
      project_id = project_id,

      output_vcf = output_vcf,
      output_vcf_index = output_vcf_index,

      target_bucket = target_bucket
  }

  output {
    String target_output_vcf = mh.target_output_vcf
    String target_output_vcf_index = mh.target_output_vcf_index
    Int target_vcf_moved=mh.target_vcf_moved
  }
}

task MoveHaplotypecallerResult {
  input {
    String genoset
    String GRID
    String? project_id

    String output_vcf
    String output_vcf_index

    String target_bucket
  }

  String new_output_vcf = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf)}"
  String new_output_vcf_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf_index)}"

  command <<<

gsutil -m ~{"-u " + project_id} mv ~{output_vcf} \
  ~{output_vcf_index} \
  ~{target_bucket}/~{genoset}/~{GRID}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_output_vcf = "~{new_output_vcf}"
    String target_output_vcf_index = "~{new_output_vcf_index}"
    Int target_vcf_moved = 1
  }
}
