version 1.0


import "../../pipelines/wdl/reprocessing/wgs/WholeGenomeReprocessing.wdl" as WholeGenomeReprocessing
import "../../verification/VerifyExomeReprocessing.wdl" as VerifyExomeReprocessing
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestWholeGenomeReprocessing {

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
      DNASeqSingleSampleReferences references
      VariantCallingScatterSettings scatter_settings
      PapiSettings papi_settings
      File? fingerprint_genotypes_file
      File? fingerprint_genotypes_index
      File wgs_coverage_interval_list
      String cloud_provider

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call WholeGenomeReprocessing.WholeGenomeReprocessing {
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
        references = references,
        scatter_settings = scatter_settings,
        papi_settings = papi_settings,
        fingerprint_genotypes_file = fingerprint_genotypes_file,
        fingerprint_genotypes_index = fingerprint_genotypes_index,
        wgs_coverage_interval_list = wgs_coverage_interval_list,
        cloud_provider = cloud_provider
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    WholeGenomeReprocessing.output_vcf_index,
                                    WholeGenomeReprocessing.output_vcf,
                                    WholeGenomeReprocessing.validate_cram_file_report,
                                    WholeGenomeReprocessing.output_cram_md5,
                                    WholeGenomeReprocessing.output_cram_index,
                                    WholeGenomeReprocessing.output_cram,
                                    WholeGenomeReprocessing.agg_quality_distribution_pdf,
                                    WholeGenomeReprocessing.agg_insert_size_histogram_pdf,
                                    WholeGenomeReprocessing.agg_gc_bias_pdf,
                                    WholeGenomeReprocessing.calculate_read_group_checksum_md5,
                                    WholeGenomeReprocessing.selfSM,
                                    WholeGenomeReprocessing.read_group_gc_bias_pdf,
                                    ],
                                    # Array[File] outputs
                                    WholeGenomeReprocessing.unsorted_read_group_quality_distribution_pdf,
                                    WholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_pdf,
                                    WholeGenomeReprocessing.unsorted_read_group_insert_size_histogram_pdf,
                                    WholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf,
                                    WholeGenomeReprocessing.unmapped_bams,
                                    WholeGenomeReprocessing.validation_report,
                                    # File? outputs
                                    select_all([WholeGenomeReprocessing.output_bqsr_reports]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    WholeGenomeReprocessing.gvcf_detail_metrics,
                                    WholeGenomeReprocessing.gvcf_summary_metrics,
                                    WholeGenomeReprocessing.duplicate_metrics,
                                    WholeGenomeReprocessing.raw_wgs_metrics,
                                    WholeGenomeReprocessing.wgs_metrics,
                                    WholeGenomeReprocessing.agg_quality_distribution_metrics,
                                    WholeGenomeReprocessing.agg_pre_adapter_summary_metrics,
                                    WholeGenomeReprocessing.agg_pre_adapter_detail_metrics,
                                    WholeGenomeReprocessing.agg_insert_size_metrics,
                                    WholeGenomeReprocessing.agg_gc_bias_summary_metrics,
                                    WholeGenomeReprocessing.agg_gc_bias_detail_metrics,
                                    WholeGenomeReprocessing.agg_bait_bias_summary_metrics,
                                    WholeGenomeReprocessing.agg_bait_bias_detail_metrics,
                                    WholeGenomeReprocessing.agg_alignment_summary_metrics,
                                    WholeGenomeReprocessing.read_group_gc_bias_summary_metrics,
                                    WholeGenomeReprocessing.read_group_gc_bias_detail_metrics,
                                    WholeGenomeReprocessing.read_group_alignment_summary_metrics,
                                    ],
                                    # Array[File] outputs
                                    WholeGenomeReprocessing.unsorted_read_group_quality_distribution_metrics,
                                    WholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_metrics,
                                    WholeGenomeReprocessing.unsorted_read_group_insert_size_metrics,
                                    WholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics,
                                    WholeGenomeReprocessing.quality_yield_metrics,
                                    # File? outputs
                                    select_all([WholeGenomeReprocessing.fingerprint_detail_metrics]),
                                    select_all([WholeGenomeReprocessing.fingerprint_summary_metrics]),
                                    select_all([WholeGenomeReprocessing.cross_check_fingerprints_metrics]),
                                    
    ])

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCram {
          input:
            input_file = WholeGenomeReprocessing.output_cram,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCrai {
          input:
            input_file = WholeGenomeReprocessing.output_cram_index,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGvcf {
          input:
            input_file = WholeGenomeReprocessing.output_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGvcfIndex {
          input:
            input_file = WholeGenomeReprocessing.output_vcf_index,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetBams {
          input:
            input_files = WholeGenomeReprocessing.unmapped_bams,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyExomeReprocessing.VerifyReprocessing as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_cram = GetCram.truth_file, 
          test_cram = GetCram.results_file,
          truth_crai = GetCrai.truth_file, 
          test_crai = GetCrai.results_file,
          truth_gvcf = GetGvcf.truth_file, 
          test_gvcf = GetGvcf.results_file,
          truth_gvcf_index = GetGvcfIndex.truth_file, 
          test_gvcf_index = GetGvcfIndex.results_file,
          truth_bams = GetBams.truth_files, 
          test_bams = GetBams.results_files,
          done = CopyToTestResults.done
      }
    }
}