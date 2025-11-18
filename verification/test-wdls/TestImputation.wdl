version 1.0


import "../../pipelines/wdl/arrays/imputation/Imputation.wdl" as Imputation
import "../../verification/VerifyImputation.wdl" as VerifyImputation
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestImputation {

    input {
      Int chunkLength = 25000000
      Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects
      File? multi_sample_vcf
      File? multi_sample_vcf_index
      Array[File]? single_sample_vcfs
      Array[File]? single_sample_vcf_indices
      Boolean perform_extra_qc_steps = false # these are optional additional extra QC steps from Amit's group that should only be
      Float? optional_qc_max_missing
      Float? optional_qc_hwe
      File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
      Array[String] contigs
      String reference_panel_path # path to the bucket where the reference panel files are stored for all contigs
      File genetic_maps_eagle
      String output_callset_name # the output callset name
      Boolean split_output_to_single_sample = false
      Int merge_ssvcf_mem_mb = 3000 # the memory allocation for MergeSingleSampleVcfs (in mb)
      Float frac_above_maf_5_percent_well_imputed_threshold = 0.9 # require fraction of maf > 0.05 sites well imputed to be greater than this to pass
      Int chunks_fail_threshold = 1 # require fewer than this many chunks to fail in order to pass
      String vcf_suffix = ".vcf.gz"
      String vcf_index_suffix = ".vcf.gz.tbi"
      String bcf_suffix = ".bcf"
      String bcf_index_suffix =  ".bcf.csi"
      String m3vcf_suffix = ".cleaned.m3vcf.gz"

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call Imputation.Imputation {
      input:
        chunkLength = chunkLength,
        chunkOverlaps = chunkOverlaps,
        multi_sample_vcf = multi_sample_vcf,
        multi_sample_vcf_index = multi_sample_vcf_index,
        single_sample_vcfs = single_sample_vcfs,
        single_sample_vcf_indices = single_sample_vcf_indices,
        perform_extra_qc_steps = perform_extra_qc_steps,
        optional_qc_max_missing = optional_qc_max_missing,
        optional_qc_hwe = optional_qc_hwe,
        ref_dict = ref_dict,
        contigs = contigs,
        reference_panel_path = reference_panel_path,
        genetic_maps_eagle = genetic_maps_eagle,
        output_callset_name = output_callset_name,
        split_output_to_single_sample = split_output_to_single_sample,
        merge_ssvcf_mem_mb = merge_ssvcf_mem_mb,
        frac_above_maf_5_percent_well_imputed_threshold = frac_above_maf_5_percent_well_imputed_threshold,
        chunks_fail_threshold = chunks_fail_threshold,
        vcf_suffix = vcf_suffix,
        vcf_index_suffix = vcf_index_suffix,
        bcf_suffix = bcf_suffix,
        bcf_index_suffix = bcf_index_suffix,
        m3vcf_suffix = m3vcf_suffix
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    Imputation.n_failed_chunks,
                                    Imputation.failed_chunks,
                                    Imputation.chunks_info,
                                    Imputation.imputed_multisample_vcf_index,
                                    Imputation.imputed_multisample_vcf,
                                    ],
                                    flatten(select_all([
                                      Imputation.imputed_single_sample_vcfs,
                                      Imputation.imputed_single_sample_vcf_indices
                                    ]))
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    Imputation.aggregated_imputation_metrics,
                                    ]
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
        call Utilities.GetValidationInputs as GetVcf {
          input:
            input_file = Imputation.imputed_multisample_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetVcfIndex {
          input:
            input_file = Imputation.imputed_multisample_vcf_index,
            results_path = results_path,
            truth_path = truth_path
        }

        if(defined(Imputation.imputed_single_sample_vcfs)){
          call Utilities.GetValidationInputs as GetImputedVCF {
            input:
              input_files = Imputation.imputed_single_sample_vcfs,
              results_path = results_path,
              truth_path = truth_path
          }
        }

      call VerifyImputation.VerifyImputation as Verify {
        input:
          split_output_to_single_sample = split_output_to_single_sample,
          output_callset_name = output_callset_name,
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_vcf = GetVcf.truth_file, 
          test_vcf = GetVcf.results_file,
          truth_vcf_index = GetVcfIndex.truth_file, 
          test_vcf_index = GetVcfIndex.results_file,
          input_multi_sample_vcf = multi_sample_vcf,
          input_multi_sample_vcf_index = multi_sample_vcf_index,
          input_single_sample_vcfs = single_sample_vcfs,
          input_single_sample_vcfs_indices = single_sample_vcf_indices,
          single_sample_test_vcf = select_first([GetImputedVCF.results_files, []]),
          single_sample_truth_vcf = select_first([GetImputedVCF.truth_files, []]),
          done = CopyToTestResults.done
      }
    }
}