version 1.0

workflow VUMCMoveCram {
  input {
    String target_bucket
    String? project_id

    String genoset
    String GRID

    String output_cram
    String output_cram_index
    String output_cram_md5
  }

  String gcs_output_dir = sub(target_bucket, "/+$", "")

  String moved_output_cram = "~{gcs_output_dir}/~{genoset}/~{GRID}/~{basename(output_cram)}"
  String moved_output_cram_index = "~{gcs_output_dir}/~{genoset}/~{GRID}/~{basename(output_cram_index)}"
  String moved_output_cram_md5 = "~{gcs_output_dir}/~{genoset}/~{GRID}/~{basename(output_cram_md5)}"

  call MoveCram as mf {
    input:
      target_bucket = gcs_output_dir,
      project_id = project_id,

      genoset = genoset,
      GRID = GRID,

      output_cram = output_cram,
      output_cram_index = output_cram_index,
      output_cram_md5 = output_cram_md5
  }

  output {
    String target_output_cram = moved_output_cram
    String target_output_cram_index = moved_output_cram_index
    String target_output_cram_md5 = moved_output_cram_md5

    Int target_cram_moved = mf.target_cram_moved
  }
}

task MoveCram {
  input {
    String target_bucket
    String? project_id

    String genoset
    String GRID

    String output_cram
    String output_cram_index
    String output_cram_md5
  }

  command <<<

set -e

gsutil -m ~{"-u " + project_id} mv \
  ~{output_cram} \
  ~{output_cram_index} \
  ~{output_cram_md5} \
  ~{target_bucket}/~{genoset}/~{GRID}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    Int target_cram_moved = 1
  }
}
