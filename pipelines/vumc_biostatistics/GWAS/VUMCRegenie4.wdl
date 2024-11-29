version 1.0

import "../../../tasks/vumc_biostatistics/WDLUtils.wdl" as WDLUtils
import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../../../tasks/vumc_biostatistics/BioUtils.wdl" as BioUtils
import "./GWASUtils.wdl" as GWASUtils

workflow VUMCRegenie4 {
  input {
    File? qc_pgen_file
    File? qc_pvar_file
    File? qc_psam_file

    File pgen_file
    File pvar_file
    File psam_file

    File phenoFile
    String phenoColList
    Boolean is_binary_traits = false

    File covarFile
    String covarColList

    String output_prefix

    String qc_option="--mac 100 --geno 0.1 --maf 0.1 --max-maf 0.9 --hwe 1e-15 --snps-only --not-chr 23-27"
    
    String step1_option="--loocv --bsize 1000"
    Int step1_block_size=1000
    
    String step2_option="--firth --approx --pThresh 0.01 --bsize 400"

    Array[String] chromosome_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

    String? billing_gcp_project_id
    String? target_gcp_folder
  }

  if(!defined(qc_pgen_file)){
    call BioUtils.PgenQCFilter {
      input:
        input_pgen = pgen_file,
        input_pvar = pvar_file,
        input_psam = psam_file,
        output_prefix = output_prefix,
        qc_option = qc_option
    }
  }
  
  File model_pgen = select_first([qc_pgen_file, PgenQCFilter.output_pgen])
  File model_pvar = select_first([qc_pvar_file, PgenQCFilter.output_pvar])
  File model_psam = select_first([qc_psam_file, PgenQCFilter.output_psam])

  call WDLUtils.count_lines as psam_count {
    input:
      input_file = model_psam
  }
  Int num_sample = psam_count.num_lines

  call WDLUtils.count_lines as pvar_count {
    input:
      input_file = model_pvar
  }
  Int num_variant = pvar_count.num_lines

  call WDLUtils.string_to_array as pheco_list {
    input:
      str = phenoColList,
      delimiter = ","
  }
  Array[String] phenotype_names = pheco_list.arr
  Int num_phenotype = length(phenotype_names)

  call WDLUtils.string_to_array as covar_list{
    input:
      str = covarColList,
      delimiter = ","
  }
  Array[String] covar_names = covar_list.arr
  Int num_covariate = length(covar_names)

  Int num_chromosome = length(chromosome_list)

  call GWASUtils.Regenie4MemoryEstimation {
    input:
      num_sample = num_sample,
      num_variant = num_variant,
      num_phenotype = num_phenotype,
      num_chromosome = num_chromosome,
      num_covariate = num_covariate,
      num_ridge_l0 = 5,
      block_size = step1_block_size
  }

  Int step1_memory_gb = Regenie4MemoryEstimation.step1_memory_gb

  call GWASUtils.Regenie4Step1FitModel as RegenieStep1FitModel {
    input:
      input_pgen = model_pgen,
      input_pvar = model_pvar,
      input_psam = model_psam,
      phenoFile = phenoFile,
      phenoColList = phenoColList,
      is_binary_traits = is_binary_traits,
      covarFile = covarFile,
      covarColList = covarColList,
      output_prefix = output_prefix,
      step1_option = step1_option,
      memory_gb = step1_memory_gb
  }

  scatter(chromosome in chromosome_list) {
    call GWASUtils.Regenie4Step2AssociationTest as RegenieStep2AssociationTest {
      input:
        pred_list_file = RegenieStep1FitModel.pred_list_file,
        pred_loco_files = RegenieStep1FitModel.pred_loco_files,
        input_pgen = pgen_file,
        input_pvar = pvar_file,
        input_psam = psam_file,
        phenoFile = phenoFile,
        phenoColList = phenoColList,
        is_binary_traits = is_binary_traits,
        covarFile = covarFile,
        covarColList = covarColList,
        output_prefix = "~{output_prefix}.chr~{chromosome}",
        step2_option = step2_option,
        chromosome = chromosome,
        memory_gb = step1_memory_gb #chromosome level memory cost would be less than step1, use step1 memory here.
    }
  }

  Array[File] regenie_chromosome_files = flatten(RegenieStep2AssociationTest.regenie_files)

  call GWASUtils.MergeRegenieChromosomeResults {
    input:
      regenie_chromosome_files = regenie_chromosome_files,
      phenotype_names = phenotype_names,
      chromosome_list = chromosome_list,
      regenie_prefix = output_prefix,
      output_prefix = output_prefix
  }

  Array[Int] indecies = range(length(phenotype_names))
  scatter(ind in indecies){
    String phenotype_name = phenotype_names[ind]
    call GWASUtils.RegeniePlots {
      input:
        regenie_file = MergeRegenieChromosomeResults.phenotype_regenie_files[ind],
        output_prefix = "~{output_prefix}.~{phenotype_name}"
    }
  }

  if(defined(target_gcp_folder)){
    String gcs_output_dir = select_first([target_gcp_folder])

    call GcpUtils.MoveOrCopyOneFile as CopyFile1 {
      input:
        source_file = RegenieStep1FitModel.pred_list_file,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = gcs_output_dir
    }
    call GcpUtils.MoveOrCopyFileArray as CopyFile2 {
      input:
        source_files = RegenieStep1FitModel.pred_loco_files,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = gcs_output_dir
    }
    scatter(output_loco_file in CopyFile2.outputFiles) {
      String pred_loco_file = output_loco_file
    }

    call GcpUtils.MoveOrCopyFileArray as CopyFile3 {
      input:
        source_files = MergeRegenieChromosomeResults.phenotype_regenie_files,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = gcs_output_dir
    }
    scatter(afile in CopyFile3.outputFiles) {
      String phenotype_regenie_file = afile
    }

    call GcpUtils.MoveOrCopyFileArray as CopyFile4 {
      input:
        source_files = RegeniePlots.qqplot_png,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = gcs_output_dir
    }
    scatter(qpng in CopyFile4.outputFiles) {
      String pheno_qqplot_png = qpng
    }

    call GcpUtils.MoveOrCopyFileArray as CopyFile5 {
      input:
        source_files = RegeniePlots.manhattan_png,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = gcs_output_dir
    }
    scatter(mpng in CopyFile5.outputFiles) {
      String pheno_manhattan_png = mpng
    }
  }

  output {
    File pred_list_file = select_first([CopyFile1.output_file, RegenieStep1FitModel.pred_list_file])
    Array[File] pred_loco_files = select_first([pred_loco_file, RegenieStep1FitModel.pred_loco_files])

    Array[File] phenotype_regenie_files = select_first([phenotype_regenie_file, MergeRegenieChromosomeResults.phenotype_regenie_files])

    Array[File] phenotype_qqplot_png = select_first([pheno_qqplot_png, RegeniePlots.qqplot_png])
    Array[File] phenotype_manhattan_png = select_first([pheno_manhattan_png, RegeniePlots.manhattan_png])
  }
}

