version 1.0

import "../../../../../tasks/wdl/ImputationBeagleQcTasks.wdl" as tasks

workflow InputQC {
    # if this changes, update the input_qc_version value in ImputationBeagle.wdl
    String pipeline_version = "1.2.2"


    input {
        Int chunkLength = 25000000
        Int chunkOverlaps = 5000000

        File multi_sample_vcf

        File ref_dict
        Array[String] contigs # list of possible contigs that will be processed. note the workflow will not error out if any of these contigs are missing
        String reference_panel_path_prefix
        String genetic_maps_path
        String output_basename

        String? pipeline_header_line
        Float? min_dr2_for_inclusion

    }

    call tasks.QcChecks {
        input:
            vcf_input = multi_sample_vcf,
            allowed_contigs = contigs,
            ref_dict = ref_dict,
    }

    output {
        Boolean passes_qc = QcChecks.passes_qc
        String qc_messages = QcChecks.qc_messages
    }
}
