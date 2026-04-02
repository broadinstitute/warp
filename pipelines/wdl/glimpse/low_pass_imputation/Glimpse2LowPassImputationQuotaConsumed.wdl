version 1.0

import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks

workflow QuotaConsumed {
    # if this changes, update the quota_consumed_version value in Glimpse2LowPassImputation.wdl
    String pipeline_version = "0.0.1"

    input {
        Array[String] contigs

        # this is the path the a directory that contains sites vcf, sites table, and reference chunks file.  should end with a "/
        String reference_panel_prefix

        # service currently does not accept VCFs as input
        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? cram_manifest
        File fasta
        File fasta_index
        String output_basename

        File ref_dict
    }

    call tasks.CountSamples {
        input:
            vcf = multi_sample_vcf
    }

    output {
        Int quota_consumed = CountSamples.nSamples
    }
}
