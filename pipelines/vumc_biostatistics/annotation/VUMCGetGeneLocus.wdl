version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../../../tasks/vumc_biostatistics/BioUtils.wdl" as BioUtils

workflow VUMCGetGeneLocus {
  input {
    String gene_symbol

    String? billing_project_id
    String? target_gcp_folder
  }

  call BioUtils.GetGeneLocus {
    input:
      gene_symbol = gene_symbol
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyOneFile as CopyFile {
      input:
        source_file = GetGeneLocus.gene_bed,
        is_move_file = false,
        project_id = billing_project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File gene_bed = select_first([CopyFile.output_file, GetGeneLocus.gene_bed])
  }
}
