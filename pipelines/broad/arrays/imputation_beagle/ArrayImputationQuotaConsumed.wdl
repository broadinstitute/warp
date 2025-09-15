version 1.0

import "../../../../tasks/broad/ImputationTasks.wdl" as tasks

workflow QuotaConsumed {
    # if this changes, update the quota_consumed_version value in ImputationBeagle.wdl
    String pipeline_version = "1.0.6"


    input {
        Int chunkLength = 25000000
        Int chunkOverlaps = 5000000

        File multi_sample_vcf

        File ref_dict
        Array[String] contigs
        String reference_panel_path_prefix
        String genetic_maps_path
        String output_basename

        String? pipeline_header_line
    }

    call tasks.CountSamples {
        input:
            vcf = multi_sample_vcf
    }

    output {
        Int quota_consumed = CountSamples.nSamples
    }
}
