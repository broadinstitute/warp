version 1.0


import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../../../tasks/vumc_biostatistics/BioUtils.wdl" as BioUtils
import "./AgdUtils.wdl" as AgdUtils

workflow VUMCPrepareGeneGenotypeWorkflow {
  input {
    String gene_symbol

    File agd_primary_grid_file
    File input_hail_mt_path_file

    Float expect_output_vcf_bgz_size_gb = 10

    String? project_id
    String? target_gcp_folder
  }

  call BioUtils.GetGeneLocus {
    input:
      gene_symbol = gene_symbol
  }

  call AgdUtils.HailMatrixExtractRegions {
    input:
      input_hail_mt_path_file = input_hail_mt_path_file,
      expect_output_vcf_bgz_size_gb = expect_output_vcf_bgz_size_gb,
      input_bed = GetGeneLocus.output_bed,
      target_prefix = gene_symbol,
      billing_project_id = select_first([project_id])
  }

  call AgdUtils.Annovar {
    input:
      input_vcf = HailMatrixExtractRegions.output_vcf,
      target_prefix = gene_symbol
  }

  call AgdUtils.PrepareGeneGenotype {
    input:
      gene_symbol = gene_symbol,
      agd_primary_grid_file = agd_primary_grid_file,
      annovar_file = Annovar.annovar_file,
      vcf_file = HailMatrixExtractRegions.output_vcf
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopySevenFiles as CopyFile {
      input:
        source_file1 = GetGeneLocus.output_bed,
        source_file2 = HailMatrixExtractRegions.output_vcf,
        source_file3 = Annovar.annovar_file,
        source_file4 = PrepareGeneGenotype.lof_genotype_file,
        source_file5 = PrepareGeneGenotype.lof_genotype_freq_file,
        source_file6 = PrepareGeneGenotype.vuc_genotype_file,
        source_file7 = PrepareGeneGenotype.vuc_genotype_freq_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File bed_file = select_first([CopyFile.output_file1, GetGeneLocus.output_bed])
    File vcf_file = select_first([CopyFile.output_file2, HailMatrixExtractRegions.output_vcf])
    File annovar_file = select_first([CopyFile.output_file3, Annovar.annovar_file])
    File lof_genotype_file = select_first([CopyFile.output_file4, PrepareGeneGenotype.lof_genotype_file])
    File lof_genotype_freq_file = select_first([CopyFile.output_file5, PrepareGeneGenotype.lof_genotype_freq_file])
    File vuc_genotype_file = select_first([CopyFile.output_file6, PrepareGeneGenotype.vuc_genotype_file])
    File vuc_genotype_freq_file = select_first([CopyFile.output_file7, PrepareGeneGenotype.vuc_genotype_freq_file])
  }
}
