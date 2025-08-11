version 1.0


import "../../pipelines/wdl/dna_seq/germline/variant_calling/VariantCalling.wdl" as VariantCalling
import "../../verification/VerifyGvcf.wdl" as VerifyGvcf
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestVariantCalling {

    input {
      Boolean run_dragen_mode_variant_calling = false
      Boolean use_spanning_event_genotyping = true
      File calling_interval_list
      File evaluation_interval_list
      Int haplotype_scatter_count
      Int break_bands_at_multiples_of
      Float? contamination
      File input_bam
      File input_bam_index
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      File? ref_str
      File dbsnp_vcf
      File dbsnp_vcf_index
      String base_file_name
      String final_vcf_base_name
      Int agg_preemptible_tries
      Boolean make_gvcf = true
      Boolean make_bamout = false
      Boolean use_gatk3_haplotype_caller = false
      Boolean skip_reblocking = false
      Boolean use_dragen_hard_filtering = false

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String cloud_provider
    }

    meta {
      allowNestedInputs: true
    }
  
    call VariantCalling.VariantCalling {
      input:
        run_dragen_mode_variant_calling = run_dragen_mode_variant_calling,
        use_spanning_event_genotyping = use_spanning_event_genotyping,
        calling_interval_list = calling_interval_list,
        evaluation_interval_list = evaluation_interval_list,
        haplotype_scatter_count = haplotype_scatter_count,
        break_bands_at_multiples_of = break_bands_at_multiples_of,
        contamination = contamination,
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_str = ref_str,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        base_file_name = base_file_name,
        final_vcf_base_name = final_vcf_base_name,
        agg_preemptible_tries = agg_preemptible_tries,
        make_gvcf = make_gvcf,
        make_bamout = make_bamout,
        use_gatk3_haplotype_caller = use_gatk3_haplotype_caller,
        skip_reblocking = skip_reblocking,
        use_dragen_hard_filtering = use_dragen_hard_filtering,
        cloud_provider = cloud_provider
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    VariantCalling.output_vcf_index,
                                    VariantCalling.output_vcf,
                                    ],
                                    # File? outputs
                                    select_all([VariantCalling.bamout_index]),
                                    select_all([VariantCalling.bamout]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    VariantCalling.vcf_detail_metrics,
                                    VariantCalling.vcf_summary_metrics,
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
        call Utilities.GetValidationInputs as GetGvcf {
          input:
            input_file = VariantCalling.output_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGvcfIndex {
          input:
            input_file = VariantCalling.output_vcf_index,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyGvcf.VerifyGvcf as Verify {
        input:
          truth_gvcf = GetGvcf.truth_file, 
          test_gvcf = GetGvcf.results_file,
          truth_gvcf_index = GetGvcfIndex.truth_file, 
          test_gvcf_index = GetGvcfIndex.results_file,
          done = CopyToTestResults.done
      }
    }
}