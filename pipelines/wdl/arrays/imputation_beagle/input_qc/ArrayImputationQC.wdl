version 1.0

import "../../../../../tasks/wdl/ImputationBeagleQcTasks.wdl" as tasks

workflow InputQC {
    # if this changes, update the input_qc_version value in ImputationBeagle.wdl
    String pipeline_version = "1.3.1"

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
