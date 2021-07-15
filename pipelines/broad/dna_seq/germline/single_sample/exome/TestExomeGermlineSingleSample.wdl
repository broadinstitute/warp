version 1.0

import "../../../../../..pipelines/broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as ExomeGermlineSingleSample
import "../../../../../../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample

workflow TestExomeGermlineSingleSample {

  input {
    PapiSettings papi_settings
    SampleAndUnmappedBams sample_and_unmapped_bams
    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index
    File target_interval_list
    File bait_interval_list
    String bait_set_name
    Boolean provide_bam_output = false

    # These values will be determined and injected into the inputs by the scala test framework
    String truth_cloud_path
    Array[String] metrics_files_to_test
    Boolean update_truth
  }

  # Run the pipeline
  call ExomeGermlineSingleSample.ExomeGermlineSingleSample {
    input:
      sample_and_unmapped_bams = sample_and_unmapped_bams,
      references = references,
      scatter_settings = scatter_settings,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings,
      target_interval_list = target_interval_list,
      bait_interval_list = bait_interval_list,
      bait_set_name = bait_set_name,
  }
  # Check if the verification should be run or skipped
  if (!update_truth) {
    # Create a list if the tst metrics files
    Array[File] test_metrics = [] #TODO: figure out how to get the metrics files that need to be tested without explicitly defining all of them

    # Gather the truth inputs for verification
    Array[String] truth_metrics =  prefix(truth_cloud_path, metrics_files_to_test)
    File truth_cram = truth_cloud_path + basename(ExomeGermlineSingleSample.output_cram)
    File truth_crai = truth_cloud_path + basename(ExomeGermlineSingleSample.output_cram_index)
    File truth_gvcf = truth_cloud_path + basename(ExomeGermlineSingleSample.output_vcf)

    # Run the verification
    call VerifyGermlineSingleSample.VerifyGermlineSingleSample {
      input:
      truth_metrics = truth_metrics,                                 # Array[File] truth_metrics
      truth_cram = truth_cram,                                       # File truth_cram
      truth_crai = truth_crai,                                       # File truth_crai
      truth_gvcf = truth_gvcf,                                       # File truth_gvcf
      test_metrics = test_metrics,                                   # Array[File] test_metric
      test_cram = ExomeGermlineSingleSample.output_cram,             # File test_cram
      test_crai = ExomeGermlineSingleSample.output_cram_index,       # File test_crai
      test_gvcf = ExomeGermlineSingleSample.output_vcf              # File test_gvcf
    }
  }
}
  # Outputs are intentionally left unddefined to automatically propagate all outputs from all tasks (wihout needing to individually define the outputs)
