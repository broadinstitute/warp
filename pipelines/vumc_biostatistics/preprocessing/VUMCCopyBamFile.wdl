version 1.0

workflow VUMCCopyBamFile {
  input {
    String source_bam
    String source_bai
    String target_bucket
  }

  call CopyBamFile {
    input:
      source_bam = source_bam,
      source_bai = source_bai,
      target_bucket = target_bucket
  }

  output {
    String target_bam = CopyBamFile.target_bam
    String target_bai = CopyBamFile.target_bai
  }
}

task CopyBamFile {
  input {
    String source_bam
    String source_bai
    String target_bucket
  }

  String bam_name = basename(source_bam)
  String bam_file = "${target_bucket}/${bam_name}"

  String bai_name = basename(source_bai)
  String bai_file = "${target_bucket}/${bai_name}"

  command <<<
  gsutil cp ~{source_bam} ~{bam_file}
  gsutil cp ~{source_bai} ~{bai_file}
>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_bam = bam_file
    String target_bai = bai_file
  }
}