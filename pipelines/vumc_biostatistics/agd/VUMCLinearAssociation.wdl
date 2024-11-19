version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./AgdUtils.wdl" as AgdUtils

workflow VUMCLinearAssociation {
  input {
    String phename
    Float phecode
    File phenotype_file

    String genotype_name
    File genotype_file

    File agd_primary_grid_file
    File demographics_file

    File pca_file
    File phecode_map_file
    File ancestry_file

    String? project_id
    String? target_gcp_folder
  }

  call AgdUtils.LinearAssociation {
    input:
      phename = phename,
      phecode = phecode,
      phenotype_file = phenotype_file,
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

  output {
    File linear_association_file = select_first([CopyFile.output_file1, LinearAssociation.linear_association_file])
    File linear_association_report = select_first([CopyFile.output_file2, LinearAssociation.linear_association_report])
  }
}
