version 1.0


import "../../pipelines/wdl/dna_seq/somatic/single_sample/ugwgs/UltimaGenomicsWholeGenomeCramOnly.wdl" as UltimaGenomicsWholeGenomeCramOnly
import "../../verification/VerifyUltimaGenomicsWholeGenomeCramOnly.wdl" as VerifyUltimaGenomicsWholeGenomeCramOnly
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestUltimaGenomicsWholeGenomeCramOnly {

    input {
      ContaminationSites contamination_sites
      AlignmentReferences alignment_references
      VcfPostProcessing vcf_post_processing
      Array[File]? input_cram_list
      Array[File]? input_bam_list
      String base_file_name
      Float rsq_threshold = 1.0
      Int reads_per_split = 20000000
      Boolean save_bam_file = false

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call UltimaGenomicsWholeGenomeCramOnly.UltimaGenomicsWholeGenomeCramOnly {
      input:
        contamination_sites = contamination_sites,
        alignment_references = alignment_references,
        vcf_post_processing = vcf_post_processing,
        input_cram_list = input_cram_list,
        input_bam_list = input_bam_list,
        base_file_name = base_file_name,
        rsq_threshold = rsq_threshold,
        reads_per_split = reads_per_split,
        save_bam_file = save_bam_file
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    UltimaGenomicsWholeGenomeCramOnly.agg_quality_distribution_pdf,
                                    UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_pdf,
                                    UltimaGenomicsWholeGenomeCramOnly.selfSM,
                                    UltimaGenomicsWholeGenomeCramOnly.output_cram_md5,
                                    UltimaGenomicsWholeGenomeCramOnly.output_cram_index,
                                    UltimaGenomicsWholeGenomeCramOnly.output_cram,
                                    ],
                                    # File? outputs
                                    select_all([UltimaGenomicsWholeGenomeCramOnly.output_bam_index]),
                                    select_all([UltimaGenomicsWholeGenomeCramOnly.output_bam]),
                                    select_all([UltimaGenomicsWholeGenomeCramOnly.agg_alignment_summary_pdf]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    UltimaGenomicsWholeGenomeCramOnly.agg_quality_distribution_metrics,
                                    UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_summary_metrics,
                                    UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_detail_metrics,
                                    UltimaGenomicsWholeGenomeCramOnly.agg_alignment_summary_metrics,
                                    UltimaGenomicsWholeGenomeCramOnly.duplicate_metrics,
                                    UltimaGenomicsWholeGenomeCramOnly.raw_wgs_metrics,
                                    UltimaGenomicsWholeGenomeCramOnly.wgs_metrics,
                                    UltimaGenomicsWholeGenomeCramOnly.quality_yield_metrics,
                                    ],
                                    
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
            input_file = UltimaGenomicsWholeGenomeCramOnly.output_cram,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCrai {
          input:
            input_file = UltimaGenomicsWholeGenomeCramOnly.output_cram_index,
            results_path = results_path,
            truth_path = truth_path
        }

        if (defined(UltimaGenomicsWholeGenomeCramOnly.output_bam) || defined(UltimaGenomicsWholeGenomeCramOnly.output_bam_index)) {
          call Utilities.ErrorWithMessage as ErrorMessageBamsShouldBeNull {
            input:
              message = "Both output_bam and output_bam_index should be null"
          }
        }

      call VerifyUltimaGenomicsWholeGenomeCramOnly.VerifyUltimaGenomicsWholeGenomeCramOnly as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_cram = GetCram.truth_file, 
          test_cram = GetCram.results_file,
          truth_crai = GetCrai.truth_file, 
          test_crai = GetCrai.results_file,
          done = CopyToTestResults.done
      }
    }
}