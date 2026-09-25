version 1.0

import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks

workflow QuotaConsumed {
    # if this changes, update the quota_consumed_version value in ImputationBeagle.wdl
    String pipeline_version = "1.1.1"

    input {
        # user provided inputs
        File multi_sample_vcf
        String output_basename
        Float? min_dr2_for_inclusion

        # service provided inputs
        Array[String] contigs # list of possible contigs that will be processed. note the workflow will not error out if any of these contigs are missing
        String genetic_maps_path
        File ref_dict
        String reference_panel_path_prefix
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
