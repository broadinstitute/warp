version 1.0


import "../../pipelines/wdl/dna_seq/germline/single_sample/ugwgs/UltimaGenomicsWholeGenomeGermline.wdl" as UltimaGenomicsWholeGenomeGermline
import "../../verification/VerifyUltimaGenomicsWholeGenomeGermline.wdl" as VerifyUltimaGenomicsWholeGenomeGermline
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestUltimaGenomicsWholeGenomeGermline {

    input {
      ContaminationSites contamination_sites
      AlignmentReferences alignment_references
      VariantCallingSettings variant_calling_settings
      VcfPostProcessing vcf_post_processing
      Array[File]? input_cram_list
      Array[File]? input_bam_list
      String base_file_name
      Float rsq_threshold = 1.0
      Boolean merge_bam_file = true
      Boolean make_haplotype_bam = false
      Int reads_per_split = 20000000
      String filtering_model_no_gt_name = "rf_model_ignore_gt_incl_hpol_runs"

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call UltimaGenomicsWholeGenomeGermline.UltimaGenomicsWholeGenomeGermline {
      input:
        contamination_sites = contamination_sites,
        alignment_references = alignment_references,
        variant_calling_settings = variant_calling_settings,
        vcf_post_processing = vcf_post_processing,
        input_cram_list = input_cram_list,
        input_bam_list = input_bam_list,
        base_file_name = base_file_name,
        rsq_threshold = rsq_threshold,
        merge_bam_file = merge_bam_file,
        make_haplotype_bam = make_haplotype_bam,
        reads_per_split = reads_per_split,
        filtering_model_no_gt_name = filtering_model_no_gt_name
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    UltimaGenomicsWholeGenomeGermline.agg_quality_distribution_pdf,
                                    UltimaGenomicsWholeGenomeGermline.agg_gc_bias_pdf,
                                    UltimaGenomicsWholeGenomeGermline.filtered_vcf_index,
                                    UltimaGenomicsWholeGenomeGermline.filtered_vcf,
                                    UltimaGenomicsWholeGenomeGermline.selfSM,
                                    UltimaGenomicsWholeGenomeGermline.output_cram_md5,
                                    UltimaGenomicsWholeGenomeGermline.output_cram_index,
                                    UltimaGenomicsWholeGenomeGermline.output_cram,
                                    UltimaGenomicsWholeGenomeGermline.output_vcf_index,
                                    UltimaGenomicsWholeGenomeGermline.output_vcf,
                                    UltimaGenomicsWholeGenomeGermline.output_gvcf_index,
                                    UltimaGenomicsWholeGenomeGermline.output_gvcf
                                    ],
                                    # File? outputs
                                    select_all([UltimaGenomicsWholeGenomeGermline.agg_alignment_summary_pdf]),
                                    select_all([UltimaGenomicsWholeGenomeGermline.haplotype_bam_index]),
                                    select_all([UltimaGenomicsWholeGenomeGermline.haplotype_bam]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    UltimaGenomicsWholeGenomeGermline.agg_quality_distribution_metrics,
                                    UltimaGenomicsWholeGenomeGermline.agg_gc_bias_summary_metrics,
                                    UltimaGenomicsWholeGenomeGermline.agg_gc_bias_detail_metrics,
                                    UltimaGenomicsWholeGenomeGermline.agg_alignment_summary_metrics,
                                    UltimaGenomicsWholeGenomeGermline.duplicate_metrics,
                                    UltimaGenomicsWholeGenomeGermline.raw_wgs_metrics,
                                    UltimaGenomicsWholeGenomeGermline.wgs_metrics,
                                    UltimaGenomicsWholeGenomeGermline.quality_yield_metrics,
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
            input_file = UltimaGenomicsWholeGenomeGermline.output_cram,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCrai {
          input:
            input_file = UltimaGenomicsWholeGenomeGermline.output_cram_index,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetVcf {
          input:
            input_file = UltimaGenomicsWholeGenomeGermline.output_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFilteredVcf {
          input:
            input_file = UltimaGenomicsWholeGenomeGermline.filtered_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFilteredVcfIndex {
          input:
            input_file = UltimaGenomicsWholeGenomeGermline.filtered_vcf_index,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGvcf {
          input:
            input_file = UltimaGenomicsWholeGenomeGermline.output_gvcf,
            results_path = results_path,
            truth_path = truth_path
        }

        call Utilities.GetValidationInputs as GetGvcfIndex {
          input:
            input_file = UltimaGenomicsWholeGenomeGermline.output_gvcf_index,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyUltimaGenomicsWholeGenomeGermline.VerifyUltimaGenomicsWholeGenomeGermline as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_cram = GetCram.truth_file, 
          test_cram = GetCram.results_file,
          truth_crai = GetCrai.truth_file, 
          test_crai = GetCrai.results_file,
          truth_vcf = GetVcf.truth_file, 
          test_vcf = GetVcf.results_file,
          truth_filtered_vcf = GetFilteredVcf.truth_file,
          truth_filtered_vcf_index = GetFilteredVcfIndex.truth_file, 
          test_filtered_vcf = GetFilteredVcf.results_file,
          test_filtered_vcf_index = GetFilteredVcfIndex.results_file,
          truth_gvcf = GetGvcf.truth_file,
          test_gvcf = GetGvcf.results_file,
          truth_gvcf_index = GetGvcfIndex.truth_file,
          test_gvcf_index = GetGvcfIndex.results_file,
          sample_name = UltimaGenomicsWholeGenomeGermline.sample_name,
          done = CopyToTestResults.done
      }
    }
}