version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifysnM3C {

    input {
        File test_mapping_summary
        File truth_mapping_summary

        Boolean? done
    }

    call VerifyTasks.CompareCompressedTextFiles as CompareMappingSummaryMetrics {
        input:
            test_zip  = test_mapping_summary,
            truth_zip = truth_mapping_summary
    }
}