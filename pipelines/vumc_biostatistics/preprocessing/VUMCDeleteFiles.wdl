version 1.0

workflow VUMCMoveFile {
  input {
    Array[String] source_files
  }

  call DeleteFiles {
    input:
      source_files = source_files,
  }

  output {
    Array[String] target_files = []
  }
}

task DeleteFiles {
  input {
    Array[String] source_files
  }
  command <<<
  gsutil -m rm -a ~{sep=" " source_files}
>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
}