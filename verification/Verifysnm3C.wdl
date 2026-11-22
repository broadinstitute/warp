version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow Verifysnm3C {

    input {
        File test_mapping_summary
        File truth_mapping_summary

        Array[File] test_name_sorted_bam_array
        Array[File] truth_name_sorted_bam_array



        Array[File] test_all_reads_3C_contacts_array
        Array[File] truth_all_reads_3C_contacts_array

        Array[File] test_cz_array
        Array[File] truth_cz_array

        # cytozip image, used to read .cz headers and decode files when a
        # byte-comparison fails
        String cytozip_docker

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

    scatter (idx in range(length(truth_cz_array))){
        call CompareCzFiles {
            input:
                test_cz  = test_cz_array[idx],
                truth_cz = truth_cz_array[idx],
                docker   = cytozip_docker
        }
    }

}

# Byte-compares a per-cell .cz file against truth. The .cz format is
# byte-reproducible for a given cytozip version, so an exact comparison is the
# right default. A cytozip upgrade changes a few header bytes (the embedded
# format version) while leaving the data untouched, which would otherwise look
# like an unexplained failure -- so on mismatch this decodes both files back to
# ALLC and reports whether the data actually changed or only the version stamp.
task CompareCzFiles {

    input {
        File test_cz
        File truth_cz
        String docker
    }

    Float file_size = size(test_cz, "GiB") + size(truth_cz, "GiB")
    Int disk_size = ceil(file_size * 10) + 20

    command <<<
        set -uo pipefail

        if cmp -s ~{test_cz} ~{truth_cz}; then
            echo "Comparison succeeded: .cz files are byte-identical."
            exit 0
        fi

        echo "Comparison failed: .cz files differ."

        field_of() {
            czip header -I "$1" 2>/dev/null | awk -v k="$2" '$1==k{print $3}'
        }
        tv=$(field_of ~{test_cz} version);  rv=$(field_of ~{truth_cz} version)
        # the .cz header records the genome fasta path it was built from, so the
        # bytes differ if the reference is localized to a different location
        tg=$(field_of ~{test_cz} message);  rg=$(field_of ~{truth_cz} message)
        echo "    cytozip format version -- test: ${tv:-unknown}  truth: ${rv:-unknown}"
        if [ "$tg" != "$rg" ]; then
            echo "    genome path in header -- test: ${tg:-unknown}  truth: ${rg:-unknown}"
        fi

        # decode both back to plain ALLC so the data can be compared
        # independently of the container's version stamp
        czip to_bgzip -I ~{test_cz}  -O test.allc.tsv.gz  --cov_col cov --allc_format > /dev/null 2>&1
        czip to_bgzip -I ~{truth_cz} -O truth.allc.tsv.gz --cov_col cov --allc_format > /dev/null 2>&1
        gzip -dc test.allc.tsv.gz  > test.allc.tsv
        gzip -dc truth.allc.tsv.gz > truth.allc.tsv

        if cmp -s test.allc.tsv truth.allc.tsv; then
            if [ "${tv:-x}" != "${rv:-y}" ]; then
                echo "    --> DATA IS IDENTICAL; only the embedded version stamp changed (${rv} -> ${tv})."
                echo "        This is an expected cytozip version bump. Re-seed truth with a"
                echo "        workflow_dispatch run using updateTruth: true."
            elif [ "$tg" != "$rg" ]; then
                echo "    --> DATA IS IDENTICAL; only the genome path recorded in the header differs."
                echo "        The reference was localized to a different location than in the truth"
                echo "        run. Harmless, but re-seed truth to silence it."
            else
                echo "    --> DATA IS IDENTICAL but the bytes differ at the same format version"
                echo "        and the same genome path. This is unexpected; the .cz encoding may"
                echo "        have changed."
            fi
        else
            echo "    --> DATA DIFFERS. This is a real regression, not a version bump."
            echo "        First differing ALLC records:"
            diff truth.allc.tsv test.allc.tsv | head -10 | sed 's/^/        /'
        fi
        exit 1
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size + " HDD"
        cpu: 1
        memory: "4 GiB"
        preemptible: 3
    }

    output {
        String done = read_string(stdout())
    }
}
