version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "Utils.wdl" as Utils

workflow VUMCMergePgenFiles {
  input {
    Array[File] pgen_files
    Array[File] pvar_files
    Array[File] psam_files

    String output_prefix

    String? project_id
    String? target_gcp_folder
  }

  call Utils.MergePgenFiles {
    input:
      pgen_files = pgen_files,
      pvar_files = pvar_files,
      psam_files = psam_files,
      output_prefix = output_prefix
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyThreeFiles as CopyFile {
      input:
        source_file1 = MergePgenFiles.output_pgen,
        source_file2 = MergePgenFiles.output_pvar,
        source_file3 = MergePgenFiles.output_psam,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_pgen = select_first([CopyFile.output_file1, MergePgenFiles.output_pgen])
    File output_pvar = select_first([CopyFile.output_file2, MergePgenFiles.output_pvar])
    File output_psam = select_first([CopyFile.output_file3, MergePgenFiles.output_psam])
  }
}
