version 1.0

task MoveOrCopyOneFile {
  input {
    String source_file

    Boolean is_move_file = false

    String? project_id
    String target_gcp_folder
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")

  String new_file = "~{gcs_output_dir}/~{basename(source_file)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{source_file} ~{gcs_output_dir}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_file = new_file
  }
}

task MoveOrCopyTwoFiles {
  input {
    String source_file1
    String source_file2

    Boolean is_move_file = false

    String? project_id
    String target_gcp_folder
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")

  String new_file1 = "~{gcs_output_dir}/~{basename(source_file1)}"
  String new_file2 = "~{gcs_output_dir}/~{basename(source_file2)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{source_file1} ~{source_file2} ~{gcs_output_dir}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_file1 = new_file1
    String output_file2 = new_file2
  }
}

task MoveOrCopyThreeFiles {
  input {
    String source_file1
    String source_file2
    String source_file3

    Boolean is_move_file = false

    String? project_id
    String target_gcp_folder
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")

  String new_file1 = "~{gcs_output_dir}/~{basename(source_file1)}"
  String new_file2 = "~{gcs_output_dir}/~{basename(source_file2)}"
  String new_file3 = "~{gcs_output_dir}/~{basename(source_file3)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{source_file1} ~{source_file2} ~{source_file3} ~{gcs_output_dir}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_file1 = new_file1
    String output_file2 = new_file2
    String output_file3 = new_file3
  }
}

task MoveOrCopyFileArray {
  input {
    Array[String] source_files

    Boolean is_move_file = false

    String? project_id
    String target_gcp_folder
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} '~{sep="' '" source_files}' ~{gcs_output_dir}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
}
