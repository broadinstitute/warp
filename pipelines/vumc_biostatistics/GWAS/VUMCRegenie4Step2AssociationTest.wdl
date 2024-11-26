version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "./GWASUtils.wdl" as GWASUtils

workflow VUMCRegenie4Step2AssociationTest {
  input {
    File pred_list_file
    Array[File] pred_loco_files

    File pgen_file
    File pvar_file
    File psam_file

    File phenoFile
    String phenoColList
    Boolean is_binary_traits = false

    File covarFile
    String covarColList

    String step2_option = "--firth --approx --pThresh 0.01 --bsize 400"

    String output_prefix

    Int? chromosome
    String? billing_gcp_project_id
    String? target_gcp_folder
  }

  call GWASUtils.Regenie4Step2AssociationTest as RegenieStep2AssociationTest {
    input:
      pred_list_file = pred_list_file,
      pred_loco_files = pred_loco_files,
      input_pgen = pgen_file,
      input_pvar = pvar_file,
      input_psam = psam_file,
      phenoFile = phenoFile,
      phenoColList = phenoColList,
      is_binary_traits = is_binary_traits,
      covarFile = covarFile,
      covarColList = covarColList,
      output_prefix = output_prefix,
      step2_option = step2_option,
      chromosome = chromosome
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyFileArray as CopyFile {
      input:
        source_files = RegenieStep2AssociationTest.regenie_files,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
    String gcs_output_dir = sub(select_first([target_gcp_folder]), "/+$", "")
    scatter(regenie_file in RegenieStep2AssociationTest.regenie_files) {
      String target_file = gcs_output_dir + "/" + basename(regenie_file)
    }
  }

  output {
    Array[File] regenie_files = select_first([target_file, RegenieStep2AssociationTest.regenie_files])
  }
}
