version 1.0

import "../../pipelines/wdl/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl" as WholeGenomeGermlineSingleSample
import "../../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestWholeGenomeGermlineSingleSample {

  input {
    SampleAndUnmappedBams sample_and_unmapped_bams
    DNASeqSingleSampleReferences references
    DragmapReference? dragmap_reference
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list

    Boolean provide_bam_output = false
    Boolean use_gatk3_haplotype_caller = true

    Boolean dragen_functional_equivalence_mode = false
    Boolean dragen_maximum_quality_mode = false

    Boolean run_dragen_mode_variant_calling = false
    Boolean use_spanning_event_genotyping = true
    Boolean unmap_contaminant_reads = true
    Boolean perform_bqsr = true
    Boolean use_bwa_mem = true
    Boolean allow_empty_ref_alt = false
    Boolean use_dragen_hard_filtering = false
    String cloud_provider

    # These values will be determined and injected into the inputs by the scala test framework
    String truth_path
    String results_path
    Boolean update_truth
  }

  meta {
    allowNestedInputs: true
  }

  # Run the pipeline
  call WholeGenomeGermlineSingleSample.WholeGenomeGermlineSingleSample {
    input:
      sample_and_unmapped_bams           = sample_and_unmapped_bams,
      references                         = references,
      dragmap_reference                  =  dragmap_reference,
      scatter_settings                   = scatter_settings,
      papi_settings                      = papi_settings,
      fingerprint_genotypes_file         = fingerprint_genotypes_file,
      fingerprint_genotypes_index        = fingerprint_genotypes_index,
      wgs_coverage_interval_list         = wgs_coverage_interval_list,
      provide_bam_output                 = provide_bam_output,
      use_gatk3_haplotype_caller         = use_gatk3_haplotype_caller,
      dragen_functional_equivalence_mode = dragen_functional_equivalence_mode,
      dragen_maximum_quality_mode        = dragen_maximum_quality_mode,
      run_dragen_mode_variant_calling    = run_dragen_mode_variant_calling,
      use_spanning_event_genotyping      = use_spanning_event_genotyping,
      unmap_contaminant_reads            = unmap_contaminant_reads,
      perform_bqsr                       = perform_bqsr,
      use_bwa_mem                        = use_bwa_mem,
      allow_empty_ref_alt                = allow_empty_ref_alt,
      use_dragen_hard_filtering          = use_dragen_hard_filtering,
      cloud_provider                     = cloud_provider
  }

  # Collect all of the pipeline outputs into a single Array[String]
  Array[String] pipeline_outputs = flatten([
                            [ # File outputs
                            WholeGenomeGermlineSingleSample.read_group_gc_bias_pdf,
                            WholeGenomeGermlineSingleSample.selfSM,
                            WholeGenomeGermlineSingleSample.calculate_read_group_checksum_md5,
                            WholeGenomeGermlineSingleSample.agg_gc_bias_pdf,
                            WholeGenomeGermlineSingleSample.agg_insert_size_histogram_pdf,
                            WholeGenomeGermlineSingleSample.agg_quality_distribution_pdf,
                            WholeGenomeGermlineSingleSample.output_cram,
                            WholeGenomeGermlineSingleSample.output_cram_index,
                            WholeGenomeGermlineSingleSample.output_cram_md5,
                            WholeGenomeGermlineSingleSample.validate_cram_file_report,
                            WholeGenomeGermlineSingleSample.output_vcf,
                            WholeGenomeGermlineSingleSample.output_vcf_index,
                            ], # Array[File] outputs
                            WholeGenomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_pdf,
                            WholeGenomeGermlineSingleSample.unsorted_read_group_insert_size_histogram_pdf,
                            WholeGenomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_pdf,
                            WholeGenomeGermlineSingleSample.unsorted_read_group_quality_distribution_pdf,
                            # File? outputs
                            select_all([WholeGenomeGermlineSingleSample.output_bqsr_reports]),
                            select_all([WholeGenomeGermlineSingleSample.output_bam]),
                            select_all([WholeGenomeGermlineSingleSample.output_bam_index])

  ])

  # Collect all of the pipeline metrics into a single Array[String]
  Array[String] pipeline_metrics = flatten([
                              [ # File outputs
                              WholeGenomeGermlineSingleSample.read_group_alignment_summary_metrics,
                              WholeGenomeGermlineSingleSample.read_group_gc_bias_detail_metrics,
                              WholeGenomeGermlineSingleSample.read_group_gc_bias_summary_metrics,
                              WholeGenomeGermlineSingleSample.agg_alignment_summary_metrics,
                              WholeGenomeGermlineSingleSample.agg_bait_bias_detail_metrics,
                              WholeGenomeGermlineSingleSample.agg_bait_bias_summary_metrics,
                              WholeGenomeGermlineSingleSample.agg_gc_bias_detail_metrics,
                              WholeGenomeGermlineSingleSample.agg_gc_bias_summary_metrics,
                              WholeGenomeGermlineSingleSample.agg_insert_size_metrics,
                              WholeGenomeGermlineSingleSample.agg_pre_adapter_detail_metrics,
                              WholeGenomeGermlineSingleSample.agg_pre_adapter_summary_metrics,
                              WholeGenomeGermlineSingleSample.agg_quality_distribution_metrics,
                              WholeGenomeGermlineSingleSample.agg_error_summary_metrics,
                              WholeGenomeGermlineSingleSample.wgs_metrics,
                              WholeGenomeGermlineSingleSample.raw_wgs_metrics,
                              WholeGenomeGermlineSingleSample.duplicate_metrics,
                              WholeGenomeGermlineSingleSample.gvcf_summary_metrics,
                              WholeGenomeGermlineSingleSample.gvcf_detail_metrics,
                              ], # Array[File] outputs
                              WholeGenomeGermlineSingleSample.quality_yield_metrics,
                              WholeGenomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_metrics,
                              WholeGenomeGermlineSingleSample.unsorted_read_group_insert_size_metrics,
                              WholeGenomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_metrics,
                              WholeGenomeGermlineSingleSample.unsorted_read_group_quality_distribution_metrics,
                              # File? outputs
                              select_all([WholeGenomeGermlineSingleSample.cross_check_fingerprints_metrics]),
                              select_all([WholeGenomeGermlineSingleSample.fingerprint_summary_metrics]),
                              select_all([WholeGenomeGermlineSingleSample.fingerprint_detail_metrics]),
  ])

  # Copy results of pipeline to test results bucket
  call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      contamination             = WholeGenomeGermlineSingleSample.contamination,
      destination_cloud_path    = results_path
  }

  # If updating truth then copy pipeline results to truth bucket
  if (update_truth){
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      contamination             = WholeGenomeGermlineSingleSample.contamination,
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
          input_file   = WholeGenomeGermlineSingleSample.output_cram,
          results_path = results_path,
          truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCrais {
      input:
        input_file   = WholeGenomeGermlineSingleSample.output_cram_index,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetGVCFs {
      input:
        input_file   = WholeGenomeGermlineSingleSample.output_vcf,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetGVCFIndexes {
      input:
        input_file    = WholeGenomeGermlineSingleSample.output_vcf_index,
        results_path  = results_path,
        truth_path    = truth_path
    }

    # done is dummy input to force copy completion before verification
    call VerifyGermlineSingleSample.VerifyGermlineSingleSample as Verify {
      input:
        truth_metrics    = GetMetricsInputs.truth_files,
        truth_cram       = GetCrams.truth_file,
        truth_crai       = GetCrais.truth_file,
        truth_gvcf       = GetGVCFs.truth_file,
        truth_gvcf_index = GetGVCFIndexes.truth_file,
        test_metrics     = GetMetricsInputs.results_files,
        test_cram        = GetCrams.results_file,
        test_crai        = GetCrais.results_file,
        test_gvcf        = GetGVCFs.results_file,
        test_gvcf_index  = GetGVCFIndexes.results_file,
        done             = CopyToTestResults.done
    }
  }

  output {
    Array[File]? metric_comparison_report_files = Verify.metric_comparison_report_files
  }

}
