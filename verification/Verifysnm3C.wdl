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


        Array[File] test_all_reads_3C_contacts_array
        Array[File] truth_all_reads_3C_contacts_array

        Array[File] test_cz_tar_array
        Array[File] truth_cz_tar_array

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

    scatter (idx in range(length(truth_unique_reads_cgn_extraction_allc_array))){
        call VerifyTasks.CompareCompressedTextFiles as Compare_unique_reads_cgn_extraction_allc_array {
            input:
                test_zip   = test_unique_reads_cgn_extraction_allc_array[idx],
                truth_zip  = truth_unique_reads_cgn_extraction_allc_array[idx]
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

    scatter (idx in range(length(truth_cz_tar_array))){
        call CompareCzFiles {
            input:
                test_cz_tar  = test_cz_tar_array[idx],
                truth_cz_tar = truth_cz_tar_array[idx],
                docker       = cytozip_docker
        }
    }

}

# Byte-compares the .cz files inside the test and truth tars. The .cz format is
# byte-reproducible for a given cytozip version, so an exact comparison is the
# right default. A cytozip upgrade changes a few header bytes (the embedded
# format version) while leaving the data untouched, which would otherwise look
# like an unexplained failure -- so on mismatch this decodes both files back to
# ALLC and reports whether the data actually changed or only the version stamp.
task CompareCzFiles {

    input {
        File test_cz_tar
        File truth_cz_tar
        String docker
    }

    Float file_size = size(test_cz_tar, "GiB") + size(truth_cz_tar, "GiB")
    Int disk_size = ceil(file_size * 10) + 20

    command <<<
        set -uo pipefail

        mkdir -p test truth
        pigz -dc ~{test_cz_tar}  | tar -xf - -C test
        pigz -dc ~{truth_cz_tar} | tar -xf - -C truth

        version_of() {
            czip header -I "$1" 2>/dev/null | awk '/^version/{print $3}'
        }

        # decode a .cz back to plain ALLC text so data can be compared
        # independently of the container's version stamp
        decode() {
            czip to_bgzip -I "$1" -O "$2.gz" --cov_col cov --allc_format >/dev/null 2>&1
            gzip -dc "$2.gz"
        }

        fail=0

        # every cell present in truth must be present in test, and match
        for truth_cz in truth/*.cz; do
            cell=$(basename "$truth_cz" .cz)
            test_cz="test/${cell}.cz"

            if [ ! -f "$test_cz" ]; then
                echo "FAIL [$cell]: present in truth but missing from test output"
                fail=1
                continue
            fi

            if cmp -s "$test_cz" "$truth_cz"; then
                continue
            fi

            echo "FAIL [$cell]: .cz files differ"
            tv=$(version_of "$test_cz")
            rv=$(version_of "$truth_cz")
            echo "    cytozip format version -- test: ${tv:-unknown}  truth: ${rv:-unknown}"

            decode "$test_cz"  "t_${cell}.allc.tsv" > "t_${cell}.allc.tsv"
            decode "$truth_cz" "r_${cell}.allc.tsv" > "r_${cell}.allc.tsv"

            if cmp -s "t_${cell}.allc.tsv" "r_${cell}.allc.tsv"; then
                if [ "${tv:-x}" != "${rv:-y}" ]; then
                    echo "    --> DATA IS IDENTICAL; only the embedded version stamp changed (${rv} -> ${tv})."
                    echo "        This is an expected cytozip version bump. Re-seed truth with a"
                    echo "        workflow_dispatch run using updateTruth: true."
                else
                    echo "    --> DATA IS IDENTICAL but the bytes differ at the same format version."
                    echo "        This is unexpected; the .cz encoding may have changed."
                fi
            else
                echo "    --> DATA DIFFERS. This is a real regression, not a version bump."
                echo "        First differing ALLC records:"
                diff "r_${cell}.allc.tsv" "t_${cell}.allc.tsv" | head -10 | sed 's/^/        /'
            fi
            fail=1
        done

        # cells in test that truth does not know about
        for test_cz in test/*.cz; do
            cell=$(basename "$test_cz" .cz)
            if [ ! -f "truth/${cell}.cz" ]; then
                echo "FAIL [$cell]: present in test output but missing from truth"
                fail=1
            fi
        done

        if [ "$fail" -eq 0 ]; then
            echo "Comparison succeeded: all .cz files are byte-identical."
        else
            exit 1
        fi
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
