version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow Verifysnm3C {

    input {
        File test_mapping_summary
        File truth_mapping_summary

        Array[File] test_name_sorted_bam_array
        Array[File] truth_name_sorted_bam_array

        Array[File] test_unique_reads_cgn_extraction_allc_array
        Array[File] truth_unique_reads_cgn_extraction_allc_array

        Array[File] test_unique_reads_cgn_extraction_allc_extract_array
        Array[File] truth_unique_reads_cgn_extraction_allc_extract_array

        Array[File] test_all_reads_3C_contacts_array
        Array[File] truth_all_reads_3C_contacts_array

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
                truth_bam  = truth_name_sorted_bam_array[idx],
                lenient_header = true
        }
   }

    scatter (idx in range(length(truth_unique_reads_cgn_extraction_allc_array))){
        call VerifyTasks.CompareCompressedTextFiles as Compare_unique_reads_cgn_extraction_allc_array {
            input:
                test_zip   = test_unique_reads_cgn_extraction_allc_array[idx],
                truth_zip  = truth_unique_reads_cgn_extraction_allc_array[idx]
        }
    }

    scatter (idx in range(length(truth_unique_reads_cgn_extraction_allc_extract_array))){
        call VerifyTasks.CompareCompressedTextFiles as Compare_unique_reads_cgn_extraction_allc_extract_array {
            input:
                test_zip   = test_unique_reads_cgn_extraction_allc_extract_array[idx],
                truth_zip  = truth_unique_reads_cgn_extraction_allc_extract_array[idx]
        }
    }

scatter (idx in range(length(truth_all_reads_3C_contacts_array))){
        call VerifyTasks.CompareCompressedTextFiles as Compare_all_reads_3C_contacts_array {
            input:
                test_zip   = test_all_reads_3C_contacts_array[idx],
                truth_zip  = truth_all_reads_3C_contacts_array[idx]
        }
    }

    output {
        Boolean done = true
    }
}
