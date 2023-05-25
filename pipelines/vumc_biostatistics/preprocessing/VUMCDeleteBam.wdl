version 1.0

workflow VUMCDeleteBam {
  input {
    String output_bam
    String output_bam_index
    String output_bam_md5
  }

  call DeleteBam as db {
    input:
      output_bam = output_bam,
      output_bam_index = output_bam_index,
      output_bam_md5 = output_bam_md5,
  }

  output {
    String target_output_bam = ""
    String target_output_bam_index = ""
    String target_output_bam_md5 = ""
  }
}

task DeleteBam {
  input {
    String output_bam
    String output_bam_index
    String output_bam_md5
  }

  command <<<

gsutil rm ~{output_bam} ~{output_bam_index} ~{output_bam_md5} 
>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
}
