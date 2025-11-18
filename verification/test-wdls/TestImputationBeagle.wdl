version 1.0


import "../../pipelines/wdl/arrays/imputation_beagle/ImputationBeagle.wdl" as ImputationBeagle
import "../../verification/VerifyImputationBeagle.wdl" as VerifyImputationBeagle
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestImputationBeagle {

    input {
      Int chunkLength = 25000000
      Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects
      Int sample_chunk_size = 1000 # this is the number of samples that will be processed in parallel in each chunked scatter
      Float min_dr2_for_inclusion = 0.0 # minimum imputation quality (DR2) for a variant to be included in the output VCF
      
      File multi_sample_vcf
      
      File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
      Array[String] contigs # list of possible contigs that will be processed. note the workflow will not error out if any of these contigs are missing
      String reference_panel_path_prefix # path + file prefix to the bucket where the reference panel files are stored for all contigs
      String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
      String output_basename # the basename for intermediate and output files

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call ImputationBeagle.ImputationBeagle {
      input:
        chunkLength = chunkLength,
        chunkOverlaps = chunkOverlaps,
        sample_chunk_size = sample_chunk_size,
        min_dr2_for_inclusion = min_dr2_for_inclusion,
        multi_sample_vcf = multi_sample_vcf,
        ref_dict = ref_dict,
        contigs = contigs,
        reference_panel_path_prefix = reference_panel_path_prefix,
        genetic_maps_path = genetic_maps_path,
        output_basename = output_basename,
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    ImputationBeagle.imputed_multi_sample_vcf,
                                    ImputationBeagle.imputed_hom_ref_sites_only_vcf
                                    ]
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    ImputationBeagle.chunks_info,
                                    ImputationBeagle.contigs_info,
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
        call Utilities.GetValidationInputs as GetImputedMultiSampleVcf {
          input:
            input_file = ImputationBeagle.imputed_multi_sample_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetImputedSitesOnlyVcf {
          input:
            input_file = ImputationBeagle.imputed_hom_ref_sites_only_vcf,
              results_path = results_path,
              truth_path = truth_path
        }


      call VerifyImputationBeagle.VerifyImputationBeagle as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          multi_sample_truth_vcf = GetImputedMultiSampleVcf.truth_file,
          multi_sample_test_vcf = GetImputedMultiSampleVcf.results_file,
          sites_only_truth_vcf = GetImputedSitesOnlyVcf.truth_file,
          sites_only_test_vcf = GetImputedSitesOnlyVcf.results_file,
          done = CopyToTestResults.done
      }
    }
}
