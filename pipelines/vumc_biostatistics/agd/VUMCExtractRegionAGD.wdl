version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../genotype/Utils.wdl" as GenotypeUtils
import "./AgdUtils.wdl" as AgdUtils

workflow VUMCExtractRegionAGD {
  input {
    Array[File] source_pgen_files
    Array[File] source_pvar_files
    Array[File] source_psam_files

    Array[String] chromosomes

    String plink2_filter_option

    File region_bed
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

    call GenotypeUtils.ExtractPgenRegions as ExtractPgenRegions {
      input:
        source_pgen = pgen_file,
        source_pvar = pvar_file,
        source_psam = ReplaceICAIdWithGrid.output_psam,
        chromosome = chromosome,
        plink2_filter_option = plink2_filter_option,
        region_bed = region_bed
    }
  }

  call GenotypeUtils.MergePgenFiles as MergePgenFiles{
    input:
      pgen_files = ExtractPgenRegions.output_pgen,
      pvar_files = ExtractPgenRegions.output_pvar,
      psam_files = ExtractPgenRegions.output_psam,
      output_prefix = target_prefix
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
