version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow Verifysnm3C {

    input {
        File test_mapping_summary
        File truth_mapping_summary

        Array[File] test_name_sorted_bam_array
        Array[File] truth_name_sorted_bam_array

        File test_unique_reads_cgn_extraction_allc_array
        File truth_unique_reads_cgn_extraction_allc_array

        Boolean? done
    }

    call VerifyTasks.CompareCompressedTextFiles as CompareMappingSummaryMetrics {
        input:
            test_zip  = test_mapping_summary,
            truth_zip = truth_mapping_summary
    }

    scatter (idx in range(length(truth_name_sorted_bam_array))){
        call VerifyTasks.CompareBams as CompareBams {
            input:
                test_bam   = test_name_sorted_bam_array[idx],
                truth_bam  = truth_name_sorted_bam_array[idx]
        }
   }

    call VerifyTasks.CompareTabix as Compare_unique_reads_cgn_extraction_allc_array {
        input:
            test_fragment_file   = test_unique_reads_cgn_extraction_allc_array,
            truth_fragment_file  = truth_unique_reads_cgn_extraction_allc_array
    }
}
