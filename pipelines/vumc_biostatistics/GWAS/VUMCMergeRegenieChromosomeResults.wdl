version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./GWASUtils.wdl" as GWASUtils

workflow VUMCMergeRegenieChromosomeResults {
  input {
    Array[File] regenie_chromosome_files
    Array[String] phenotype_names
    Array[Int] chromosome_list

    String regenie_prefix

    String output_prefix

    String? billing_gcp_project_id
    String? target_gcp_folder
  }

  call GWASUtils.MergeRegenieChromosomeResults {
    input:
      regenie_chromosome_files = regenie_chromosome_files,
      phenotype_names = phenotype_names,
      chromosome_list = chromosome_list,
      regenie_prefix = regenie_prefix,
      output_prefix = output_prefix
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyFileArray as CopyFile {
      input:
        source_files = MergeRegenieChromosomeResults.phenotype_regenie_files,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
    scatter(regenie_file in CopyFile.outputFiles) {
      String target_file = regenie_file
    }
  }

  output {
    Array[File] phenotype_regenie_files = select_first([target_file, MergeRegenieChromosomeResults.phenotype_regenie_files])
  }
}
