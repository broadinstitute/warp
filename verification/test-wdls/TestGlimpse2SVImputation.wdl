version 1.0

import "../../pipelines/wdl/glimpse/sv_imputation/Glimpse2SVImputation.wdl" as Glimpse2SVImputation
import "../../verification/VerifyGlimpse2SVImputation.wdl" as VerifyGlimpse2SVImputation
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestGlimpse2SVImputation {

    input {
        # inputs for Preprocessing wdl
        File? input_gvcfs_fofn
        File? input_gvcf_idxs_fofn
        File? sample_ids_file          # order of sample ids must match that of gVCFs

        Array[File]? input_gvcfs
        Array[File]? input_gvcf_idxs
        Array[String]? sample_ids

        String output_prefix

        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        String? extract_bubble_likelihoods_extra_args

        Array[String] paste_regions

        # inputs for Batch wdl
        Array[String] chromosomes
        File genetic_maps_tsv
        File ref_dict
        File chunked_panel_json

        String? extra_phase_args

        # override for cpu used for glimpse phase task. Mostly used to set to 1 for determinism in testing, defaults to 4
        Int? glimpse_phase_cpu_override

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json

        String? glimpse2_docker

        # These values will be determined and injected into the inputs by the scala test framework
        String truth_path
        String results_path
        Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }

    call Glimpse2SVImputation.Glimpse2SVImputation {
      input:
        input_gvcfs_fofn = input_gvcfs_fofn,
        input_gvcf_idxs_fofn = input_gvcf_idxs_fofn,
        sample_ids_file = sample_ids_file,
        input_gvcfs = input_gvcfs,
        input_gvcf_idxs = input_gvcf_idxs,
        sample_ids = sample_ids,
        output_prefix = output_prefix,
        preprocess_panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
        preprocess_panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
        extract_bubble_likelihoods_extra_args = extract_bubble_likelihoods_extra_args,
        paste_regions = paste_regions,
        chromosomes = chromosomes,
        genetic_maps_tsv = genetic_maps_tsv,
        chunked_panel_json = chunked_panel_json,
        ref_dict = ref_dict,
        extra_phase_args = extra_phase_args,
        glimpse_phase_cpu_override = glimpse_phase_cpu_override,
        pop_glimpse2_panel_resources_json = pop_glimpse2_panel_resources_json,
        glimpse2_docker = glimpse2_docker,
    }


    # Collect all pipeline outputs into a single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    Glimpse2SVImputation.glimpse2_popped_posteriors_vcf,
                                    Glimpse2SVImputation.glimpse2_popped_posteriors_vcf_idx
    ])

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = pipeline_outputs,
        destination_cloud_path    = results_path
    }

    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input:
          files_to_copy             = pipeline_outputs,
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetPoppedPosteriorsVcfs {
          input:
            input_files = Glimpse2SVImputation.glimpse2_popped_posteriors_vcf,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyGlimpse2SVImputation.VerifyGlimpse2SVImputation as Verify {
        input:
          truth_popped_posteriors_vcf = GetPoppedPosteriorsVcfs.truth_files,
          test_popped_posteriors_vcf = GetPoppedPosteriorsVcfs.results_files,
          done = CopyToTestResults.done
      }
    }
}
