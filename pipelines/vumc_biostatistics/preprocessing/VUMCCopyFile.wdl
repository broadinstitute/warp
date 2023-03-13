version 1.0

workflow VUMCCopyFile {
  input {
    String source_file
    String? project_id
    String target_bucket
  }

  call CopyFile {
    input:
      source_file = source_file,
      project_id = project_id,
      target_bucket = target_bucket
  }

  output {
    String target_file = CopyFile.target_file
    Array[String] target_file_array = [CopyFile.target_file]
  }
}

task CopyFile {
  input {
    String source_file
    String? project_id
    String target_bucket
  }

  String target_url = "${target_bucket}/${basename(source_file)}"

  command <<<
  gsutil ~{"-u " + project_id} cp ~{source_file} ~{target_url}
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