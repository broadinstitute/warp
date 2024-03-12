version 1.0

import "./Utils.wdl" as Utils

workflow VUMCVcfMoveFiles {
  input {
    String input_vcf
    String input_vcf_index

    String? project_id
    String target_gcp_folder
  }

  call Utils.MoveOrCopyTwoFiles as MoveFiles {
    input:
      source_file1 = input_vcf,
      source_file2 = input_vcf_index,
      is_move_file = true,
      project_id = project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    File output_vcf = MoveFiles.output_file1
    File output_vcf_index = MoveFiles.output_file2
  }
}

