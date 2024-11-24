version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./AgdUtils.wdl" as AgdUtils

workflow VUMCPhewas {
  input {
    File phecode_list_file
    Int min_occurance = 2

    String genotype_name
    File genotype_file

    File agd_primary_grid_file
    File pca_file

    File demographics_file
    File ancestry_file

    File phecode_data_file
    File phecode_map_file

    String? project_id
    String? target_gcp_folder
  }

  Array[String] phecodes = read_lines(phecode_list_file)

  scatter (phecode in phecodes){
    call AgdUtils.PreparePhenotype {
      input:
        phename = "~{phecode}",
        phecode = phecode,
        agd_primary_grid_file = agd_primary_grid_file,
        phecode_data_file = phecode_data_file,
        phecode_map_file = phecode_map_file,
        min_occurance = min_occurance
    }

    call AgdUtils.LinearAssociation {
      input:
        phename = "~{phecode}",
        phecode = phecode,
        phenotype_file = PreparePhenotype.phenotype_file,
        agd_primary_grid_file = agd_primary_grid_file,
        demographics_file = demographics_file,
        genotype_name = genotype_name,
        genotype_file = genotype_file,
        pca_file = pca_file,
        phecode_map_file = phecode_map_file,
        ancestry_file = ancestry_file
    }

    if(defined(target_gcp_folder)){
      call GcpUtils.MoveOrCopyTwoFiles as CopyFile {
        input:
          source_file1 = LinearAssociation.linear_association_file,
          source_file2 = LinearAssociation.linear_association_report,
          is_move_file = false,
          project_id = project_id,
          target_gcp_folder = select_first([target_gcp_folder])
      }
    }
    File linear_association_file = select_first([CopyFile.output_file1, LinearAssociation.linear_association_file])
    File linear_association_report = select_first([CopyFile.output_file2, LinearAssociation.linear_association_report])
  }

  call AgdUtils.LinearAssociationSummary {
    input:
      phecode_map_file = phecode_map_file,
      phecodes = phecodes,
      linear_association_files = LinearAssociation.linear_association_file,
      genotype_name = genotype_name
  }

  output {
    Array[File] linear_association_files = linear_association_file
    Array[File] linear_association_reports = linear_association_report
    File linear_association_summary_file = LinearAssociationSummary.linear_association_summary_file
  }
}
