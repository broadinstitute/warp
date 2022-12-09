version 1.0

workflow VUMCMoveFile {
  input {
    String source_file
    String target_bucket
  }

  call MoveFile {
    input:
      source_file = source_file,
      target_bucket = target_bucket
  }

  output {
    String target_file = MoveFile.target_file
  }
}

task MoveFile {
  input {
    String source_file
    String target_bucket
  }

  String target_url = "${target_bucket}/${basename(source_file)}"

  command <<<
  gsutil mv ~{source_file} ~{target_url}
>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_file = "~{target_url}"
  }
}