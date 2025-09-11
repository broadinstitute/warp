version 1.0

import "../../../../tasks/broad/ImputationBeagleQcTasks.wdl" as tasks

workflow InputQC {
    String pipeline_version = "1.0.3" # if this changes, update the input_qc_version value in ImputationBeagle.wdl as well

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

    call tasks.QcChecks {
        input:
            vcf_input = multi_sample_vcf,
            ref_dict = ref_dict,
    }

    output {
        Boolean passes_qc = QcChecks.passes_qc
        String qc_messages = QcChecks.qc_messages
    }
}
