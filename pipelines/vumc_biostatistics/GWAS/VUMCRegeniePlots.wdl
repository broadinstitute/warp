version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./GWASUtils.wdl" as GWASUtils

workflow VUMCRegeniePlots {
  input {
    File regenie_file
    String output_prefix

    String? billing_gcp_project_id
    String? target_gcp_folder
  }

  call GWASUtils.RegeniePlots {
    input:
      regenie_file = regenie_file,
      output_prefix = output_prefix
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyTwoFiles as CopyFile {
      input:
        source_file1 = RegeniePlots.qqplot_png,
        source_file2 = RegeniePlots.manhattan_png,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File qqplot_png = select_first([CopyFile.output_file1, RegeniePlots.qqplot_png])
    File manhattan_png = select_first([CopyFile.output_file2, RegeniePlots.manhattan_png])
  }
}
