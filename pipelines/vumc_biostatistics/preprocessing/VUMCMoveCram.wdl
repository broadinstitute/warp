version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils

workflow VUMCMoveCram {
  input {
    String target_bucket
    String? project_id

    String genoset
    String GRID

    String input_cram
    String input_cram_index
    String input_cram_md5
  }

  String gcs_output_dir = sub(target_bucket, "/+$", "")
  String target_gcp_folder = "~{gcs_output_dir}/~{genoset}/~{GRID}"

  call GcpUtils.MoveOrCopyThreeFiles as CopyFile {
    input:
      source_file1 = input_cram,
      source_file2 = input_cram_index,
      source_file3 = input_cram_md5,

      is_move_file = false,

      project_id = project_id,
      target_gcp_folder = target_gcp_folder,
  }

  output {
    String output_cram = CopyFile.output_file1
    String output_cram_index = CopyFile.output_file2
    String output_cram_md5 = CopyFile.output_file3
    Int output_cram_moved = 1
  }
}
