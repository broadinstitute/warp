version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./AgdUtils.wdl" as AgdUtils

workflow VUMCPreparePhenotype {
  input {
    Float phecode

    File agd_primary_grid_file
    File phecode_data_file
    File phecode_map_file
    Int min_occurance = 2

    String? project_id
    String? target_gcp_folder
  }

  call AgdUtils.PreparePhenotype {
    input:
      phecode = phecode,
      agd_primary_grid_file = agd_primary_grid_file,
      phecode_data_file = phecode_data_file,
      phecode_map_file = phecode_map_file,
      min_occurance = min_occurance
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyTwoFiles as CopyFile {
      input:
        source_file1 = PreparePhenotype.phenotype_file,
        source_file2 = PreparePhenotype.phenotype_report,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File phenotype_file = select_first([CopyFile.output_file1, PreparePhenotype.phenotype_file])
    File phenotype_report = select_first([CopyFile.output_file2, PreparePhenotype.phenotype_report])
  }
}
