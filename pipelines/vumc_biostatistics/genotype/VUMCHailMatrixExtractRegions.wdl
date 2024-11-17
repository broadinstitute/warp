version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../agd/AgdUtils.wdl" as AgdUtils

workflow VUMCHailMatrixExtractRegions {
  input {
    File input_hail_mt_path_file
    File input_bed

    Float expect_output_vcf_bgz_size_gb

    String target_prefix

    String billing_project_id
    String? target_gcp_folder
  }

  call AgdUtils.HailMatrixExtractRegions {
    input:
      input_hail_mt_path_file = input_hail_mt_path_file,
      expect_output_vcf_bgz_size_gb = expect_output_vcf_bgz_size_gb,
      input_bed = input_bed,
      target_prefix = target_prefix,
      billing_project_id = billing_project_id
  }

  call AgdUtils.Annovar {
    input:
      input_vcf = HailMatrixExtractRegions.output_vcf,
      target_prefix = target_prefix
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyTwoFiles as CopyFile {
      input:
        source_file1 = HailMatrixExtractRegions.output_vcf,
        source_file2 = Annovar.output_annovar_file,
        is_move_file = false,
        project_id = billing_project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_vcf = select_first([CopyFile.output_file1, HailMatrixExtractRegions.output_vcf])
    File output_annovar_file = select_first([CopyFile.output_file2, Annovar.output_annovar_file])
  }
}
