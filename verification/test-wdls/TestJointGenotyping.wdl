version 1.0


import "../../pipelines/wdl/dna_seq/germline/joint_genotyping/JointGenotyping.wdl" as JointGenotyping
import "../../verification/VerifyJointGenotyping.wdl" as VerifyJointGenotyping
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestJointGenotyping {

    input {
      File unpadded_intervals_file
      String callset_name
      File sample_name_map
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      File dbsnp_vcf
      File dbsnp_vcf_index
      Int small_disk
      Int medium_disk
      Int large_disk
      Int huge_disk
      Array[String]? snp_recalibration_tranche_values
      Array[String] snp_recalibration_annotation_values
      Array[String]? indel_recalibration_tranche_values
      Array[String]? indel_recalibration_annotation_values
      File haplotype_database
      File eval_interval_list
      File hapmap_resource_vcf
      File hapmap_resource_vcf_index
      File omni_resource_vcf
      File omni_resource_vcf_index
      File one_thousand_genomes_resource_vcf
      File one_thousand_genomes_resource_vcf_index
      File mills_resource_vcf
      File mills_resource_vcf_index
      File axiomPoly_resource_vcf
      File axiomPoly_resource_vcf_index
      File dbsnp_resource_vcf = dbsnp_vcf
      File dbsnp_resource_vcf_index = dbsnp_vcf_index
      Float excess_het_threshold = 54.69
      Float? vqsr_snp_filter_level
      Float? vqsr_indel_filter_level
      File? targets_interval_list
      Int? snp_vqsr_downsampleFactor
      Int? top_level_scatter_count
      Boolean? gather_vcfs
      Boolean? run_vets
      Int snps_variant_recalibration_threshold = 500000
      Boolean rename_gvcf_samples = true
      Float unbounded_scatter_count_scale_factor = 0.15
      Int gnarly_scatter_count = 10
      Boolean use_gnarly_genotyper = false
      Boolean use_allele_specific_annotations = true
      Boolean cross_check_fingerprints = true
      Boolean scatter_cross_check_fingerprints = false

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call JointGenotyping.JointGenotyping {
      input:
        unpadded_intervals_file = unpadded_intervals_file,
        callset_name = callset_name,
        sample_name_map = sample_name_map,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        small_disk = small_disk,
        medium_disk = medium_disk,
        large_disk = large_disk,
        huge_disk = huge_disk,
        snp_recalibration_tranche_values = snp_recalibration_tranche_values,
        snp_recalibration_annotation_values = snp_recalibration_annotation_values,
        indel_recalibration_tranche_values = indel_recalibration_tranche_values,
        indel_recalibration_annotation_values = indel_recalibration_annotation_values,
        haplotype_database = haplotype_database,
        eval_interval_list = eval_interval_list,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        mills_resource_vcf = mills_resource_vcf,
        mills_resource_vcf_index = mills_resource_vcf_index,
        axiomPoly_resource_vcf = axiomPoly_resource_vcf,
        axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        excess_het_threshold = excess_het_threshold,
        vqsr_snp_filter_level = vqsr_snp_filter_level,
        vqsr_indel_filter_level = vqsr_indel_filter_level,
        snp_vqsr_downsampleFactor = snp_vqsr_downsampleFactor,
        targets_interval_list = targets_interval_list,
        top_level_scatter_count = top_level_scatter_count,
        gather_vcfs = gather_vcfs,
        run_vets = run_vets,
        snps_variant_recalibration_threshold = snps_variant_recalibration_threshold,
        rename_gvcf_samples = rename_gvcf_samples,
        unbounded_scatter_count_scale_factor = unbounded_scatter_count_scale_factor,
        gnarly_scatter_count = gnarly_scatter_count,
        use_gnarly_genotyper = use_gnarly_genotyper,
        use_allele_specific_annotations = use_allele_specific_annotations,
        cross_check_fingerprints = cross_check_fingerprints,
        scatter_cross_check_fingerprints = scatter_cross_check_fingerprints
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    # Array[File] outputs
                                    JointGenotyping.output_intervals,
                                    JointGenotyping.output_vcf_indices,
                                    JointGenotyping.output_vcfs,
                                    # File outputs
                                    select_all([
                                    JointGenotyping.crosscheck_fingerprint_check,
                                    ])
    ])


    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    JointGenotyping.summary_metrics_file,
                                    JointGenotyping.detail_metrics_file,
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
        call Utilities.GetValidationInputs as GetVcfs {
          input:
            input_files = JointGenotyping.output_vcfs,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetVcfIndexes {
          input:
            input_files = JointGenotyping.output_vcf_indices,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetIntervals {
          input:
            input_files = JointGenotyping.output_intervals,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }

        if (cross_check_fingerprints){
        call Utilities.GetValidationInputs as GetFingerprint {
            input:
              input_file = JointGenotyping.crosscheck_fingerprint_check,
              results_path = results_path,
              truth_path = truth_path
          }
        }

      call VerifyJointGenotyping.VerifyJointGenotyping as Verify {
        input:
          truth_vcfs = GetVcfs.truth_files, 
          test_vcfs = GetVcfs.results_files,
          truth_vcf_indexes = GetVcfIndexes.truth_files, 
          test_vcf_indexes = GetVcfIndexes.results_files,
          truth_intervals = GetIntervals.truth_files, 
          test_intervals = GetIntervals.results_files,
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_fingerprint = GetFingerprint.truth_file,
          test_fingerprint = GetFingerprint.results_file,
          done = CopyToTestResults.done
      }
    }
}