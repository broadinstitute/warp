version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./AgdUtils.wdl" as AgdUtils

workflow VUMCPrepareGeneGenotype {
  input {
    String gene_symbol
    File agd_primary_grid_file
    File annovar_file
    File vcf_file

    Int loss_of_function_only=1

    String? project_id
    String? target_gcp_folder
  }

  call AgdUtils.PrepareGeneGenotype {
    input:
      gene_symbol = gene_symbol,
      agd_primary_grid_file = agd_primary_grid_file,
      annovar_file = annovar_file,
      vcf_file = vcf_file,
      loss_of_function_only = loss_of_function_only
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyOneFile as CopyFile {
      input:
        source_file = PrepareGeneGenotype.genotype_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File genotype_file = select_first([CopyFile.output_file, PrepareGeneGenotype.genotype_file])
  }
}
