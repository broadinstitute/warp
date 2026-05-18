version 1.0

import "../../pipelines/wdl/glimpse/low_pass_imputation/Glimpse2LowPassImputation.wdl" as Glimpse2LowPassImputation
import "../../verification/VerifyGlimpse2LowPassImputation.wdl" as VerifyGlimpse2LowPassImputation
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestGlimpse2LowPassImputation {

    input {
        String reference_panel_prefix

        Array[String] contigs

        File? cram_manifest
        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        
        Int sample_batch_size = 1000

        Int? glimpse_phase_cpu_override


        File fasta
        File fasta_index
        String output_basename

        File ref_dict

        # These values will be determined and injected into the inputs by the scala test framework
        String truth_path
        String results_path
        Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call Glimpse2LowPassImputation.Glimpse2LowPassImputation {
      input:
        contigs = contigs,
        reference_panel_prefix = reference_panel_prefix,
        cram_manifest = cram_manifest,
        crams = crams,
        cram_indices = cram_indices,
        sample_ids = sample_ids,
        fasta = fasta,
        fasta_index = fasta_index,
        ref_dict = ref_dict,
        sample_batch_size = sample_batch_size,
        glimpse_phase_cpu_override = glimpse_phase_cpu_override,
        output_basename = output_basename,
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    Glimpse2LowPassImputation.imputed_vcf,
                                    Glimpse2LowPassImputation.imputed_hom_ref_sites_only_vcf
                                    ]
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    Glimpse2LowPassImputation.qc_metrics
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
            input_file = Glimpse2LowPassImputation.imputed_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetImputedSitesOnlyVcf {
          input:
            input_file = Glimpse2LowPassImputation.imputed_hom_ref_sites_only_vcf,
            results_path = results_path,
            truth_path = truth_path
        }


      call VerifyGlimpse2LowPassImputation.VerifyGlimpse2LowPassImputation as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          multi_sample_truth_vcf = GetImputedMultiSampleVcf.truth_file,
          multi_sample_test_vcf = GetImputedMultiSampleVcf.results_file,
          hom_ref_truth_vcf = GetImputedSitesOnlyVcf.truth_file,
          hom_ref_test_vcf = GetImputedSitesOnlyVcf.results_file,
          done = CopyToTestResults.done
      }
    }
}
