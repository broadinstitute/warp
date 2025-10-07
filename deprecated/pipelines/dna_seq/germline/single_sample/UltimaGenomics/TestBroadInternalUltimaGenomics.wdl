version 1.0


import "../../pipelines/wdl/internal/dna_seq/germline/single_sample/UltimaGenomics/BroadInternalUltimaGenomics.wdl" as BroadInternalUltimaGenomics
import "../../verification/VerifyUltimaGenomicsWholeGenomeGermline.wdl" as VerifyUltimaGenomicsWholeGenomeGermline
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestBroadInternalUltimaGenomics {

    input {
      String sample_lsid
      String output_basename
      File haplotype_database_file = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"

      String? tdr_dataset_uuid

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
      String collab_sample_id_run_id

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String vault_token_path
      String vault_token_path_arrays
      String google_account_vault_path
      String environment
    }

    meta {
      allowNestedInputs: true
    }
  
    call BroadInternalUltimaGenomics.BroadInternalUltimaGenomics {
      input:
        environment = environment,
        sample_lsid = sample_lsid,
        output_basename = output_basename,
        haplotype_database_file = haplotype_database_file,
        tdr_dataset_uuid = tdr_dataset_uuid,
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
        filtering_model_no_gt_name = filtering_model_no_gt_name,
        collab_sample_id_run_id = collab_sample_id_run_id,
        vault_token_path = vault_token_path_arrays,
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    BroadInternalUltimaGenomics.agg_quality_distribution_pdf,
                                    BroadInternalUltimaGenomics.agg_gc_bias_pdf,
                                    BroadInternalUltimaGenomics.filtered_vcf_index,
                                    BroadInternalUltimaGenomics.filtered_vcf,
                                    BroadInternalUltimaGenomics.selfSM,
                                    BroadInternalUltimaGenomics.output_cram_md5,
                                    BroadInternalUltimaGenomics.output_cram_index,
                                    BroadInternalUltimaGenomics.output_cram,
                                    BroadInternalUltimaGenomics.output_vcf_index,
                                    BroadInternalUltimaGenomics.output_vcf,
                                    BroadInternalUltimaGenomics.output_gvcf_index,
                                    BroadInternalUltimaGenomics.output_gvcf,
                                    ],
                                    # File? outputs
                                    select_all([BroadInternalUltimaGenomics.agg_alignment_summary_pdf]),
                                    select_all([BroadInternalUltimaGenomics.haplotype_bam_index]),
                                    select_all([BroadInternalUltimaGenomics.haplotype_bam]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    BroadInternalUltimaGenomics.agg_quality_distribution_metrics,
                                    BroadInternalUltimaGenomics.agg_gc_bias_summary_metrics,
                                    BroadInternalUltimaGenomics.agg_gc_bias_detail_metrics,
                                    BroadInternalUltimaGenomics.agg_alignment_summary_metrics,
                                    BroadInternalUltimaGenomics.duplicate_metrics,
                                    BroadInternalUltimaGenomics.raw_wgs_metrics,
                                    BroadInternalUltimaGenomics.wgs_metrics,
                                    BroadInternalUltimaGenomics.quality_yield_metrics,
                                    ],
                                    # File? outputs
                                    select_all([BroadInternalUltimaGenomics.picard_fingerprint_detail_metrics]),
                                    select_all([BroadInternalUltimaGenomics.picard_fingerprint_summary_metrics]),
                                    
    ])

    # Copy results of pipeline to test results bucket
    call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        vault_token_path          = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
          vault_token_path          = vault_token_path,
          google_account_vault_path = google_account_vault_path,
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
            input_file = BroadInternalUltimaGenomics.output_cram,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCrai {
          input:
            input_file = BroadInternalUltimaGenomics.output_cram_index,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetVcf {
          input:
            input_file = BroadInternalUltimaGenomics.output_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFilteredVcf {
          input:
            input_file = BroadInternalUltimaGenomics.filtered_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFilteredVcfIndex {
          input:
            input_file = BroadInternalUltimaGenomics.filtered_vcf_index,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGvcf {
          input:
            input_file = BroadInternalUltimaGenomics.output_gvcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGvcfIndex {
            input:
                input_file = BroadInternalUltimaGenomics.output_gvcf_index,
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
          sample_name = BroadInternalUltimaGenomics.sample_name,
          done = CopyToTestResults.done
      }
    }
}