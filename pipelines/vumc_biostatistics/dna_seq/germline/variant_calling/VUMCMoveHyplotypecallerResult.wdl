version 1.0

workflow VUMCMoveHyplotypecallerResult {
  input {
    String genoset
    String GRID

    String output_vcf
    String output_vcf_index

    String target_bucket
  }

  call MoveHyplotypecallerResult as mh {
    input:
      genoset = genoset,
      GRID = GRID,

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

task MoveHyplotypecallerResult {
  input {
    String genoset
    String GRID

    String output_vcf
    String output_vcf_index

    String target_bucket
  }

  String new_output_vcf = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf)}"
  String new_output_vcf_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf_index)}"

  command <<<
  gsutil mv ~{output_vcf} ~{new_output_vcf}
  gsutil mv ~{output_vcf_index} ~{new_output_vcf_index}
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
