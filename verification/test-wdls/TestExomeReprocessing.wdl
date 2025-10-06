version 1.0

import "../../pipelines/wdl/reprocessing/exome/ExomeReprocessing.wdl" as ExomeReprocessing
import "../../verification/VerifyExomeReprocessing.wdl" as VerifyExomeReprocessing
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy
import "../../structs/dna_seq/DNASeqStructs.wdl"


workflow TestExomeReprocessing {

  input {
    File? input_cram
    File? input_bam
    File? output_map
    String sample_name

    String base_file_name

    String final_gvcf_base_name

    String unmapped_bam_suffix

    File? cram_ref_fasta
    File? cram_ref_fasta_index
    String cloud_provider

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index
    File target_interval_list
    File bait_interval_list
    String bait_set_name

    # These values will be determined and injected into the inputs by the scala test framework
    String truth_path
    String results_path
    Boolean update_truth
  }


  meta {
    allowNestedInputs: true
  }

  # Run the pipeline
  call ExomeReprocessing.ExomeReprocessing {
    input:
    input_cram = input_cram,
    input_bam = input_bam,
    output_map = output_map,
    sample_name = sample_name,
    base_file_name = base_file_name,
    final_gvcf_base_name = final_gvcf_base_name,
    unmapped_bam_suffix = unmapped_bam_suffix,
    cram_ref_fasta = cram_ref_fasta,
    cram_ref_fasta_index = cram_ref_fasta_index,
    fingerprint_genotypes_file = fingerprint_genotypes_file,
    fingerprint_genotypes_index = fingerprint_genotypes_index,
    target_interval_list = target_interval_list,
    bait_interval_list = bait_interval_list,
    bait_set_name = bait_set_name,
    references = references,
    scatter_settings = scatter_settings,
    papi_settings = papi_settings,
    cloud_provider = cloud_provider
  }

  # Collect all of the pipeline outputs into a single Array[String]]
  Array[String] pipeline_outputs = flatten([
                            [ # File outputs
                            ExomeReprocessing.selfSM,
                            ExomeReprocessing.calculate_read_group_checksum_md5,
                            ExomeReprocessing.output_cram,
                            ExomeReprocessing.output_cram_index,
                            ExomeReprocessing.output_cram_md5,
                            ExomeReprocessing.validate_cram_file_report,
                            ExomeReprocessing.output_vcf,
                            ExomeReprocessing.output_vcf_index,
                            ExomeReprocessing.agg_quality_distribution_pdf,
                            ExomeReprocessing.agg_insert_size_histogram_pdf,
                            ], # Array[File] outputs
                            ExomeReprocessing.validation_report,
                            ExomeReprocessing.unmapped_bams,
                            ExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf,
                            ExomeReprocessing.unsorted_read_group_insert_size_histogram_pdf,
                            ExomeReprocessing.unsorted_read_group_quality_by_cycle_pdf,
                            ExomeReprocessing.unsorted_read_group_quality_distribution_pdf,
                            # File? outputs
                            select_all([ExomeReprocessing.output_bqsr_reports]),
  ])

  # Collect all of the pipeline metrics into a single Array[String]
  Array[String] pipeline_metrics = flatten([
                              [ # File outputs
                              ExomeReprocessing.read_group_alignment_summary_metrics,
                              ExomeReprocessing.agg_alignment_summary_metrics,
                              ExomeReprocessing.agg_bait_bias_detail_metrics,
                              ExomeReprocessing.agg_bait_bias_summary_metrics,
                              ExomeReprocessing.agg_insert_size_metrics,
                              ExomeReprocessing.agg_pre_adapter_detail_metrics,
                              ExomeReprocessing.agg_pre_adapter_summary_metrics,
                              ExomeReprocessing.agg_quality_distribution_metrics,
                              ExomeReprocessing.agg_error_summary_metrics,
                              ExomeReprocessing.duplicate_metrics,
                              ExomeReprocessing.gvcf_summary_metrics,
                              ExomeReprocessing.gvcf_detail_metrics,
                              ExomeReprocessing.hybrid_selection_metrics,
                              ], # Array[File] outputs
                              ExomeReprocessing.quality_yield_metrics,
                              ExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics,
                              ExomeReprocessing.unsorted_read_group_insert_size_metrics,
                              ExomeReprocessing.unsorted_read_group_quality_by_cycle_metrics,
                              ExomeReprocessing.unsorted_read_group_quality_distribution_metrics,
                              # File? outputs
                              select_all([ExomeReprocessing.cross_check_fingerprints_metrics]),
                              select_all([ExomeReprocessing.fingerprint_summary_metrics]),
                              select_all([ExomeReprocessing.fingerprint_detail_metrics]),
  ])

  # Copy results of pipeline to test results bucket
  call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      destination_cloud_path    = results_path
  }

  # If updating truth then copy pipeline results to truth bucket
  if (update_truth){
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      destination_cloud_path    = truth_path
    }
  }

  # If not updating truth then we need to collect all input for the validation WDL
  # This is achieved by passing each desired file/array[files] to GetValidationInputs
  if (!update_truth){
    call Utilities.GetValidationInputs as GetMetricsInputs {
      input:
        input_files  = pipeline_metrics,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCrams {
      input:
        input_file   = ExomeReprocessing.output_cram,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetBams {
      input:
        input_files  = ExomeReprocessing.unmapped_bams,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCrais {
      input:
        input_file   = ExomeReprocessing.output_cram_index,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetGVCFs {
      input:
        input_file   = ExomeReprocessing.output_vcf,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetGVCFIndexes {
      input:
        input_file    = ExomeReprocessing.output_vcf_index,
        results_path  = results_path,
        truth_path    = truth_path
    }


    # done is dummy input to force copy completion before verification
    call VerifyExomeReprocessing.VerifyReprocessing as Verify {
      input:
        test_bams = GetBams.results_files,
        truth_bams = GetBams.truth_files,
        truth_metrics = GetMetricsInputs.truth_files,
        test_metrics = GetMetricsInputs.results_files,
        truth_cram = GetCrams.truth_file,
        truth_crai = GetCrais.truth_file,
        test_cram = GetCrams.results_file,
        test_crai = GetCrais.results_file,
        truth_gvcf = GetGVCFs.truth_file,
        truth_gvcf_index = GetGVCFIndexes.truth_file,
        test_gvcf = GetGVCFs.results_file,
        test_gvcf_index = GetGVCFIndexes.results_file,
        done = CopyToTestResults.done
    }
  }


}