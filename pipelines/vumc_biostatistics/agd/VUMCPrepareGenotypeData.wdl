version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../genotype/Utils.wdl" as GenotypeUtils
import "./AgdUtils.wdl" as AgdUtils

workflow VUMCPrepareGenotypeData {
  input {
    Array[File] source_pgen_files
    Array[File] source_pvar_files
    Array[File] source_psam_files

    Array[String] chromosomes

    String plink2_filter_option

    File grid_file
    String target_prefix

    File id_map_file

    String? project_id
    String? target_gcp_folder
  }

  scatter (idx in range(length(chromosomes))) {
    String chromosome = chromosomes[idx]
    File pgen_file = source_pgen_files[idx]
    File pvar_file = source_pvar_files[idx]
    File psam_file = source_psam_files[idx]
    String replaced_sample_name = "~{chromosome}.psam"

    call AgdUtils.ReplaceICAIdWithGrid as ReplaceICAIdWithGrid {
      input:
        input_psam = psam_file,
        id_map_file = id_map_file,
        output_psam = replaced_sample_name
    }

    String grid_sample_name = "~{chromosome}.grid.psam"
    call  AgdUtils.CreateCohortPsam as CreateCohortPsam {
      input:
        input_psam = ReplaceICAIdWithGrid.output_psam,
        grid_file = grid_file,
        output_psam = grid_sample_name
    }

    call GenotypeUtils.ExtractPgenSamples as ExtractPgenSamples {
      input:
        source_pgen = pgen_file,
        source_pvar = pvar_file,
        source_psam = ReplaceICAIdWithGrid.output_psam,
        chromosome = chromosome,
        plink2_filter_option = plink2_filter_option,
        extract_sample = CreateCohortPsam.output_psam
    }
  }

  call GenotypeUtils.MergePgenFiles as MergePgenFiles{
    input:
      pgen_files = ExtractPgenSamples.output_pgen_file,
      pvar_files = ExtractPgenSamples.output_pvar_file,
      psam_files = ExtractPgenSamples.output_psam_file,
      target_prefix = target_prefix
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyThreeFiles as CopyFile {
      input:
        source_file1 = MergePgenFiles.output_pgen_file,
        source_file2 = MergePgenFiles.output_pvar_file,
        source_file3 = MergePgenFiles.output_psam_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_pgen_file = select_first([CopyFile.output_file1, MergePgenFiles.output_pgen_file])
    File output_pvar_file = select_first([CopyFile.output_file2, MergePgenFiles.output_pvar_file])
    File output_psam_file = select_first([CopyFile.output_file3, MergePgenFiles.output_psam_file])
  }
}
