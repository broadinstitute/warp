version 1.0

import "../../../../../../pipelines/broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as ExomeGermlineSingleSample
import "../../../../../../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample
import "../../../../../../tasks/broad/Utilities.wdl" as Utilities


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
    String? truth_path
    String? results_path
    Boolean? use_timestamp
    Boolean? update_truth
    String? timestamp
    String? cromwell_env
    #Array[String] metrics_files_to_test
    #Boolean update_truth
  }

  meta {
    allowNestedInputs: true
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

  call Utilities.CopyWorfklowOutputsByPath as CopyToTestResults {
    input:
      output_file_path = ExomeGermlineSingleSample.output_vcf,
      copy_bucket_path = results_path,
      workflow_name = "ExomeGermlineSingleSample",
      cromwell_env = cromwell_env,
  }



}
  # Outputs are intentionally left undefined to automatically propagate all outputs from all tasks (without needing to individually define the outputs)
