version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./VUMCPrepareGeneGenotypeWorkflow.wdl" as doGeneGenotype
import "./VUMCPhewas.wdl" as doPhewas

workflow VUMCGenePhewas {
  input {
    String gene_symbol

    File agd_primary_grid_file
    File input_hail_mt_path_file

    Float expect_output_vcf_bgz_size_gb = 10

    File phecode_list_file
    Int min_occurance = 2

    File pca_file

    File demographics_file
    File ancestry_file

    File phecode_data_file
    File phecode_map_file

    String billing_gcp_project_id
    String? target_gcp_folder
  }

  call doGeneGenotype.VUMCPrepareGeneGenotypeWorkflow as GeneGenotype {
    input:
      gene_symbol = gene_symbol,
      agd_primary_grid_file = agd_primary_grid_file,
      input_hail_mt_path_file = input_hail_mt_path_file,
      expect_output_vcf_bgz_size_gb = expect_output_vcf_bgz_size_gb,
      project_id = billing_gcp_project_id,
      target_gcp_folder = target_gcp_folder
  }

  call doPhewas.VUMCPhewas as Phewas {
    input:
      phecode_list_file = phecode_list_file,
      min_occurance = min_occurance,
      genotype_name = gene_symbol + ".lof",
      genotype_file = GeneGenotype.lof_genotype_file,
      agd_primary_grid_file = agd_primary_grid_file,
      pca_file = pca_file,
      demographics_file = demographics_file,
      ancestry_file = ancestry_file,
      phecode_data_file = phecode_data_file,
      phecode_map_file = phecode_map_file,
      project_id = billing_gcp_project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    File gene_bed = GeneGenotype.gene_bed
    File vcf_file = GeneGenotype.vcf_file
    File annovar_file = GeneGenotype.annovar_file
    String lof_genotype_name = GeneGenotype.lof_genotype_name
    File lof_genotype_file = GeneGenotype.lof_genotype_file
    File lof_genotype_freq_file = GeneGenotype.lof_genotype_freq_file
    String vus_genotype_name = GeneGenotype.vus_genotype_name
    File vus_genotype_file = GeneGenotype.vus_genotype_file
    File vus_genotype_freq_file = GeneGenotype.vus_genotype_freq_file

    Array[File] phenotype_files = Phewas.phenotype_files
    Array[File] phenotype_reports = Phewas.phenotype_reports
    Array[File] linear_association_files = Phewas.linear_association_files
    Array[File] linear_association_reports = Phewas.linear_association_reports
    File linear_association_summary_file = Phewas.linear_association_summary_file
  }
}
