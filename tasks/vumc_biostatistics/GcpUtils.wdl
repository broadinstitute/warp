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

task MoveOrCopyFourFiles {
  input {
    String source_file1
    String source_file2
    String source_file3
    String source_file4

    Boolean is_move_file = false

    String? project_id
    String target_gcp_folder
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")

  String new_file1 = "~{gcs_output_dir}/~{basename(source_file1)}"
  String new_file2 = "~{gcs_output_dir}/~{basename(source_file2)}"
  String new_file3 = "~{gcs_output_dir}/~{basename(source_file3)}"
  String new_file4 = "~{gcs_output_dir}/~{basename(source_file4)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{source_file1} ~{source_file2} ~{source_file3} ~{source_file4} ~{gcs_output_dir}/

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
    String output_file4 = new_file4
  }
}


task MoveOrCopySevenFiles {
  input {
    String source_file1
    String source_file2
    String source_file3
    String source_file4
    String source_file5
    String source_file6
    String source_file7

    Boolean is_move_file = false

    String? project_id
    String target_gcp_folder
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")

  String new_file1 = "~{gcs_output_dir}/~{basename(source_file1)}"
  String new_file2 = "~{gcs_output_dir}/~{basename(source_file2)}"
  String new_file3 = "~{gcs_output_dir}/~{basename(source_file3)}"
  String new_file4 = "~{gcs_output_dir}/~{basename(source_file4)}"
  String new_file5 = "~{gcs_output_dir}/~{basename(source_file5)}"
  String new_file6 = "~{gcs_output_dir}/~{basename(source_file6)}"
  String new_file7 = "~{gcs_output_dir}/~{basename(source_file7)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{source_file1} ~{source_file2} ~{source_file3} ~{source_file4} ~{source_file5} ~{source_file6} ~{source_file7} ~{gcs_output_dir}/

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
    String output_file4 = new_file4
    String output_file5 = new_file5
    String output_file6 = new_file6
    String output_file7 = new_file7
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

  scatter(file in source_files) {
    String new_file = "~{gcs_output_dir}/~{basename(file)}"
  }

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

  output {
    Array[String] output_files = new_file
  }
}
