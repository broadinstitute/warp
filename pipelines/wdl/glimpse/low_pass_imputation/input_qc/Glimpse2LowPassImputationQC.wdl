version 1.0

workflow InputQC {
    # if this changes, update the input_qc_version value in Glimpse2LowPassImputation.wdl
    String pipeline_version = "1.1.0"

    input {
        # service expects only cram_manifest even though main wdl can alternatively take input arrays
        File cram_manifest
        String output_basename
        Float? info_filter_for_inclusion

        Array[String] contigs
        # this is the path to a directory that contains sites vcf, sites table, and reference chunks file. should end with a "/"
        String reference_panel_prefix
        File fasta
        File fasta_index
        File ref_dict

        # used for warp tests only (which use inputs in an RP bucket). service does not support RP buckets and will not provide this input.
        String? billing_project_for_rp

        # optional additional header line to add to the output VCF
        String? pipeline_header_line
    }

    call ConvertCramManifestToInputArrays {
        input:
            cram_manifest = cram_manifest
    }

    if (ConvertCramManifestToInputArrays.passes_qc) {
        call ValidateCramsAndIndicesAndSampleIds {
            input:
                crams = ConvertCramManifestToInputArrays.crams,
                cram_indices = ConvertCramManifestToInputArrays.cram_indices,
                billing_project_for_rp = billing_project_for_rp
        }
    }

    # only check cram contents if the previous QC checks passed
    if (ConvertCramManifestToInputArrays.passes_qc && ValidateCramsAndIndicesAndSampleIds.passes_qc) {
        call ValidateCramContents {
            input:
                crams = ConvertCramManifestToInputArrays.crams,
                contigs = contigs,
                ref_dict = ref_dict,
                billing_project_for_rp = billing_project_for_rp
        }
    }

    output {
        Boolean passes_qc = select_first([ValidateCramContents.passes_qc, ValidateCramsAndIndicesAndSampleIds.passes_qc, ConvertCramManifestToInputArrays.passes_qc])
        String qc_messages = select_first([ValidateCramContents.qc_messages, ValidateCramsAndIndicesAndSampleIds.qc_messages, ConvertCramManifestToInputArrays.qc_messages])
    }
}


task ConvertCramManifestToInputArrays {
    input {
        File cram_manifest
    }

    command <<<
        cat <<EOF > script.py
        import pandas as pd

        qc_messages = []

        qc_messages_filename = "qc_messages.txt"
        passes_qc_filename = "passes_qc.txt"
        crams_filename = "crams.txt"
        cram_indices_filename = "cram_indices.txt"

        def write_column(column_data, filename):
            """Write column to file, with each value stripped of leading/trailing whitespace."""
            filtered = column_data.fillna('').astype(str).str.strip()
            with open(filename, 'w') as f:
                for value in filtered:
                    f.write(f"{value}\n")

        # Read the manifest
        try:
            df = pd.read_csv("~{cram_manifest}", sep='\t')

            # Check for required columns
            required_cols = ['cram_path', 'cram_index_path']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                with open(qc_messages_filename, 'w') as qc_file:
                    qc_file.write(f"Missing required column header(s) in the CRAM manifest: {', '.join(missing_cols)}.")
                with open(passes_qc_filename, 'w') as f:
                    f.write("false")

                # Create empty output files
                open(crams_filename, 'w').close()
                open(cram_indices_filename, 'w').close()
            else:
                # Write to output files, stripping leading/trailing whitespace from each value
                write_column(df['cram_path'], crams_filename)
                write_column(df['cram_index_path'], cram_indices_filename)

                # Write QC results
                with open(qc_messages_filename, 'w') as f:
                    f.write('\n'.join(qc_messages) if qc_messages else '')

                with open(passes_qc_filename, 'w') as f:
                    f.write("true" if not qc_messages else "false")

        except Exception as e:
            with open(qc_messages_filename, 'w') as qc_file:
                qc_file.write(f"Error reading CRAM manifest: {str(e)}.")
            with open(passes_qc_filename, 'w') as f:
                f.write("false")

            # Create empty output files
            open(crams_filename, 'w').close()
            open(cram_indices_filename, 'w').close()

        EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "4 GiB"
        preemptible: 3
        maxRetries: 2
        noAddress: true
    }

    output {
        Array[String] crams = read_lines("crams.txt")
        Array[String] cram_indices = read_lines("cram_indices.txt")
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}


task ValidateCramsAndIndicesAndSampleIds {
    input {
        Array[String] crams
        Array[String] cram_indices

        Int max_cram_file_size_gb = 10
        String? billing_project_for_rp # if set, will use this to check file sizes for requester pays buckets. if not set and input is in a RP bucket, and check will fail
    }

    String billing_project = select_first([billing_project_for_rp, ""])

    command <<<
        pip install google-cloud-storage

        cat <<'EOF' > script.py
        from google.cloud import storage

        qc_messages = []

        # Parse WDL arrays from space-separated strings
        parsed_crams = """~{sep=' ' crams}""".split()
        parsed_cram_indices = """~{sep=' ' cram_indices}""".split()

        # remove empty strings
        crams = [c for c in parsed_crams if c]
        cram_indices = [c for c in parsed_cram_indices if c]

        num_crams = len(crams)
        num_cram_indices = len(cram_indices)

        MAX_ITEMS_IN_ERROR_MESSAGES = 5

        def create_error_message_with_item_list(base_error_message: str, items_list: list) -> str:
            """Helper function to create error messages that include a list of items, but truncates the list if it's too long.
            """
            items_list_to_show = items_list[:MAX_ITEMS_IN_ERROR_MESSAGES]
            exceeded_limit_message = f"; first {MAX_ITEMS_IN_ERROR_MESSAGES} are" if len(items_list) > MAX_ITEMS_IN_ERROR_MESSAGES else ""
            return f"{base_error_message}{exceeded_limit_message}: {', '.join(items_list_to_show)}."

        def pluralize(number: int, subject: str) -> str:
            """Helper function to return a properly pluralized phrase based on the number provided, e.g. '1 CRAM file' or '2 CRAM files'."""
            return f"{number} {subject}" if number == 1 else f"{number} {subject}s"

        # Validate that the number of CRAMs and CRAIs match
        if num_crams != num_cram_indices:
            qc_messages.append(f"Found different numbers of CRAMs ({num_crams}) and CRAM index files ({num_cram_indices}).")
        else:
            print(f"Number of CRAMs and CRAM index files match: found {num_crams} of each.")

        # Ensure all crams end with .cram and all cram indices end with .crai
        crams_with_wrong_extension = [c for c in crams if not c.endswith('.cram')]
        if crams_with_wrong_extension:
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(crams_with_wrong_extension), 'CRAM file')} that do not have a .cram extension", crams_with_wrong_extension))
        else:
            print("All CRAM files have the correct .cram extension.")

        cram_indices_with_wrong_extension = [c for c in cram_indices if not c.endswith('.cram.crai')]
        if cram_indices_with_wrong_extension:
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(cram_indices_with_wrong_extension), 'CRAM index file')} without a .cram.crai extension", cram_indices_with_wrong_extension))
        else:
            print("All CRAM index files have the correct .cram.crai extension.")

        # Validate that cram paths are unique
        unique_crams = set(crams)
        if len(unique_crams) != num_crams:
            duplicates = [c for c in unique_crams if crams.count(c) > 1]
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(duplicates), 'set')} of duplicate CRAM paths", duplicates))
        else:
            print("CRAM paths are unique.")

        # Validate that each cram-crai pair has matching basenames
        mismatched_basename_pairs = []
        for cram, crai in zip(crams, cram_indices):
            cram_basename = cram.split('/')[-1].removesuffix('.cram')
            crai_basename = crai.split('/')[-1].removesuffix('.cram.crai')
            if cram_basename != crai_basename:
                mismatched_basename_pairs.append(f"{cram} and {crai}")
        if mismatched_basename_pairs:
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(mismatched_basename_pairs), 'CRAM-CRAI pair')} with mismatched basenames", mismatched_basename_pairs))
        else:
            print("All CRAM-CRAI pairs have matching basenames.")

        # Ensure that all CRAM files are less than the maximum file size allowed
        max_cram_file_size_gb = ~{max_cram_file_size_gb}
        billing_project = "~{billing_project}"
        print(f"Using billing project '{billing_project}' to check file sizes for requester pays buckets." if billing_project else "No billing project provided for requester pays buckets; file size checks may fail for files in requester pays buckets.")

        crams_exceeding_max_size = []
        files_with_access_issues = []
        files_with_invalid_gcs_format = []

        client = storage.Client()
        for cram in crams:
            try:
                if cram.startswith('gs://'):
                    # Parse GCS URI
                    bucket_name, blob_name = cram[5:].split('/', 1)

                    # Get bucket and set user_project for requester pays
                    if billing_project:
                        bucket = client.bucket(bucket_name, user_project=billing_project)
                    else:
                        bucket = client.bucket(bucket_name)

                    blob = bucket.blob(blob_name)

                    # Reload to get metadata
                    blob.reload(client=client)

                    file_size_bytes = int(blob.size)
                    file_size_gb = file_size_bytes // (1024 ** 3)
                    print(f" - File size for {cram}: {file_size_gb} GB")

                    if file_size_gb > max_cram_file_size_gb:
                        crams_exceeding_max_size.append(f"{cram} ({file_size_gb}GB)")
                else:
                    files_with_invalid_gcs_format.append(cram)
            except Exception as e:
                files_with_access_issues.append(cram)
                print(f"ERROR DETAILS for file size check for {cram}: {str(e)}")

        for crai in cram_indices:
            try:
                if crai.startswith('gs://'):
                    # Parse GCS URI
                    bucket_name, blob_name = crai[5:].split('/', 1)

                    # Get bucket and set user_project for requester pays
                    if billing_project:
                        bucket = client.bucket(bucket_name, user_project=billing_project)
                    else:
                        bucket = client.bucket(bucket_name)

                    blob = bucket.blob(blob_name)

                    # Reload to ensure the file exists and is accessible
                    blob.reload(client=client)
                else:
                    files_with_invalid_gcs_format.append(crai)
            except Exception as e:
                files_with_access_issues.append(crai)
                print(f"ERROR DETAILS for file size check for {crai}: {str(e)}")

        if files_with_invalid_gcs_format:
            qc_messages.append(create_error_message_with_item_list(
                f"Found {pluralize(len(files_with_invalid_gcs_format), 'file')} with invalid GCS format (must start with 'gs://')",
                files_with_invalid_gcs_format))

        if files_with_access_issues:
            qc_messages.append(create_error_message_with_item_list(
                f"Found {pluralize(len(files_with_access_issues), 'file')} that could not be accessed (may be due to non-existent files, lack of permissions, or requester pays bucket)",
                files_with_access_issues))

        if crams_exceeding_max_size:
            qc_messages.append(create_error_message_with_item_list(
                f"Found {pluralize(len(crams_exceeding_max_size), 'CRAM file')} exceeding the maximum allowed file size of {max_cram_file_size_gb}GB",
                crams_exceeding_max_size))
        else:
            print(f"All CRAM files are within the maximum allowed file size of {max_cram_file_size_gb}GB.")

        # Write output files
        with open("qc_messages.txt", 'w') as f:
            f.write('\n'.join(qc_messages) if qc_messages else '')

        with open("passes_qc.txt", 'w') as f:
            f.write("true" if not qc_messages else "false")

        EOF
        python3 script.py
    >>>

    runtime {
        docker:  "us.gcr.io/broad-gatk/gatk:4.6.1.0" # has python 3.10
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "4 GiB"
        maxRetries: 2
    }

    output {
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}

task ValidateCramContents {
    input {
        Array[String] crams
        Array[String] contigs
        File ref_dict
        String? billing_project_for_rp
        Int cpu = 4
    }

    String billing_project = select_first([billing_project_for_rp, ""])
    String ref_dict_basename = basename(ref_dict)

    command <<<
        # set up auth for accessing files using samtools
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # configure billing project to use for requester pays buckets, if billing project provided
        if [ -n "~{billing_project}" ]; then
            echo "Using billing project '~{billing_project}' for requester pays buckets."
            export GCS_REQUESTER_PAYS_PROJECT=~{billing_project}
        fi

        touch qc_messages.txt

        contigs=(~{sep=' ' contigs})
        ref_dict="~{ref_dict}"
        ref_dict_basename="~{ref_dict_basename}"

        declare -A ref_md5sums
        expected_count=${#contigs[@]}
        found_count=0

        while read -r line; do
            if [[ $line == @SQ* ]]; then
                # chrom is in the 2nd column of the ref dict in format SN:<chromName>
                chrom=$(echo "$line" | cut -f 2 | cut -d ":" -f 2)
                # md5sum is in the 4th column of the ref dict in format M5:<md5sum>
                md5=$(echo "$line" | cut -f 4 | cut -d ":" -f 2)

                if [[ " ${contigs[@]} " =~ " ${chrom} " ]]; then
                    ref_md5sums["$chrom"]="$md5"
                    found_count=$((found_count + 1))

                    if [[ $found_count -eq $expected_count ]]; then
                        echo "Found all ${found_count} expected contigs in reference dictionary."
                        break
                    fi
                fi
            fi
        done < $ref_dict

        echo "found relevant contigs with these md5sums in ref dict ${ref_dict}:"
        for chrom in "${!ref_md5sums[@]}"; do
            echo "  $chrom: ${ref_md5sums[$chrom]}"
        done

        MAX_ITEMS_IN_ERROR_MESSAGES=5
        cpu_count=~{cpu}

        # Split the full CRAM list into $cpu_count round-robin chunks, one per worker, so we can
        # validate all CRAMs in parallel using the CPUs available to this VM instead of checking
        # them one at a time. `r/N` distributes lines round-robin across chunks (rather than
        # contiguous line ranges), gracefully creating empty chunk files when there are fewer
        # CRAMs than workers; `-d`/`--additional-suffix` keep the chunks/chunk_* glob below working.
        printf '%s\n' ~{sep=' ' crams} > all_crams.txt
        mkdir -p chunks results
        split -n "r/${cpu_count}" -d --additional-suffix=.txt all_crams.txt chunks/chunk_

        # Validates every CRAM listed in $1 (one path per line), writing this worker's list of
        # problem CRAMs for each check to results/<worker_id>_<check>.txt. Stops early once this
        # worker's own chunk has already accumulated more than MAX_ITEMS_IN_ERROR_MESSAGES issues,
        # since the final aggregated message (built after all workers finish) is truncated to that
        # many examples anyway.
        check_cram_chunk() {
            local chunk_file="$1"
            local worker_id="$2"

            local worker_crams_with_bad_or_missing_md5sums=()
            local worker_crams_with_multiple_samples=()
            local worker_cram_sample_ids=()

            while IFS= read -r cram; do
                [ -z "$cram" ] && continue
                echo "[worker $worker_id] Validating CRAM file: $cram"
                header=$(samtools view -H "$cram")

                # check that cram is single-sample, and record its sample ID(s) so we can check for
                # sample IDs duplicated across CRAMs once all workers finish. -P (Perl regex) is
                # needed for two reasons here: (1) without it, `\t` inside a `[^\t]` bracket
                # expression isn't recognized as an actual tab, so the match would run to the end of
                # the line instead of stopping at the next field, swallowing subsequent @RG fields
                # (e.g. LB, PL) into the sample ID; (2) it enables \K, which discards everything
                # matched before it from the output, so -o prints just the ID after "SM:" without a
                # separate sed strip.
                mapfile -t sample_ids_in_cram < <(echo "$header" | grep '^@RG' | grep -oP 'SM:\K[^\t]*' | sort -u)
                n_samples=${#sample_ids_in_cram[@]}
                if [ "$n_samples" -ne 1 ]; then
                    echo "[worker $worker_id] CRAM file $cram contains data for $n_samples samples; expected exactly 1."
                    worker_crams_with_multiple_samples+=("$cram")
                fi
                # record each (sample ID, CRAM) pair so that once all workers finish, we can report
                # not just which sample IDs are duplicated across CRAMs but which CRAMs they came from
                for sample_id in "${sample_ids_in_cram[@]}"; do
                    worker_cram_sample_ids+=("$sample_id"$'\t'"$cram")
                done

                cram_ok=true
                for chrom in "${!ref_md5sums[@]}"; do
                    expected_md5=${ref_md5sums[$chrom]}
                    if ! echo "$header" | grep -q "SN:$chrom.*M5:$expected_md5"; then
                        echo "[worker $worker_id] CRAM file $cram is missing expected reference alignment MD5 for contig $chrom or it does not match the expected value."
                        worker_crams_with_bad_or_missing_md5sums+=("$cram")
                        cram_ok=false
                        break # no need to check other contigs for this cram if one is already missing or has a bad md5sum
                    fi
                done
                if [ "$cram_ok" = true ]; then
                    echo "[worker $worker_id] CRAM file $cram contains expected reference alignment MD5sums for all expected contigs"
                fi

                # stop early once this worker's own chunk already has enough issues to fill a truncated message
                total_issue_count=$(( ${#worker_crams_with_bad_or_missing_md5sums[@]} + ${#worker_crams_with_multiple_samples[@]} ))
                if [ "$total_issue_count" -gt "$MAX_ITEMS_IN_ERROR_MESSAGES" ]; then
                    echo "[worker $worker_id] found more than $MAX_ITEMS_IN_ERROR_MESSAGES CRAM files with issues in this chunk; skipping the rest of this worker's chunk"
                    break
                fi
            done < "$chunk_file"

            # Write each result list to its own file, one path per line (or leave the file empty).
            # The length is checked before expanding "${arr[@]}", since expanding a zero-element
            # array directly is unsafe under `set -u` on some older bash versions.
            if [ ${#worker_crams_with_bad_or_missing_md5sums[@]} -gt 0 ]; then
                printf '%s\n' "${worker_crams_with_bad_or_missing_md5sums[@]}" > "results/${worker_id}_bad_md5sum.txt"
            else
                : > "results/${worker_id}_bad_md5sum.txt"
            fi
            if [ ${#worker_crams_with_multiple_samples[@]} -gt 0 ]; then
                printf '%s\n' "${worker_crams_with_multiple_samples[@]}" > "results/${worker_id}_multi_sample.txt"
            else
                : > "results/${worker_id}_multi_sample.txt"
            fi
            if [ ${#worker_cram_sample_ids[@]} -gt 0 ]; then
                printf '%s\n' "${worker_cram_sample_ids[@]}" > "results/${worker_id}_sample_ids.txt"
            else
                : > "results/${worker_id}_sample_ids.txt"
            fi
        }

        worker_id=0
        for chunk_file in chunks/chunk_*; do
            check_cram_chunk "$chunk_file" "$worker_id" &
            worker_id=$((worker_id + 1))
        done
        wait

        # Merge every worker's partial results back into single lists before applying the final,
        # truncated aggregate message (same truncation behavior as before, just applied once at the
        # end instead of while looping through CRAMs one at a time).
        mapfile -t crams_with_bad_or_missing_md5sums < <(cat results/*_bad_md5sum.txt 2>/dev/null)
        mapfile -t crams_with_multiple_samples < <(cat results/*_multi_sample.txt 2>/dev/null)
        mapfile -t all_cram_sample_id_pairs < <(cat results/*_sample_ids.txt 2>/dev/null)

        # Group the (sample ID, CRAM) pairs by sample ID, and report only sample IDs associated with
        # more than one distinct CRAM, formatted as "sample_id (cram_1, cram_2)". `sort -u` first
        # collapses exact-duplicate pairs (e.g. the same CRAM's ID appearing twice due to multiple
        # @RG lines) so the count reflects distinct CRAMs, not repeated mentions of the same one.
        mapfile -t duplicate_sample_id_descriptions < <(
            printf '%s\n' "${all_cram_sample_id_pairs[@]}" | sort -u | awk -F'\t' '
                { crams[$1] = (crams[$1] == "" ? $2 : crams[$1] ", " $2); count[$1]++ }
                END { for (id in count) if (count[id] > 1) print id " (" crams[id] ")" }
            ' | sort
        )

        # if crams_with_bad_or_missing_md5sums is not empty, write an error message to qc_messages.txt
        n_bad_crams=${#crams_with_bad_or_missing_md5sums[@]}
        if [ $n_bad_crams -ne 0 ]; then
            {
                # Show only first N items if list is too long
                if [ $n_bad_crams -gt $MAX_ITEMS_IN_ERROR_MESSAGES ]; then
                    first_part_of_message="Found more than $MAX_ITEMS_IN_ERROR_MESSAGES CRAM files not aligned to the expected reference ($ref_dict_basename)"
                    second_part_of_message="; first $MAX_ITEMS_IN_ERROR_MESSAGES are:"
                    joined=$(IFS=","; echo "${crams_with_bad_or_missing_md5sums[*]:0:$MAX_ITEMS_IN_ERROR_MESSAGES}")
                    list_to_show="${joined//,/, }" # Replaces every ',' with ', '
                    echo "$first_part_of_message$second_part_of_message $list_to_show"
                else
                    if [ $n_bad_crams -eq 1 ]; then
                        pluralized=""
                    else
                        pluralized="s"
                    fi
                    first_part_of_message="Found $n_bad_crams CRAM file$pluralized not aligned to the expected reference ($ref_dict_basename)"
                    joined=$(IFS=","; echo "${crams_with_bad_or_missing_md5sums[*]}")
                    list_to_show="${joined//,/, }" # Replaces every ',' with ', '
                    echo "$first_part_of_message: $list_to_show"
                fi
            } >> qc_messages.txt
        else
            echo "All CRAM files contain the expected reference alignment MD5sums for the expected contigs."
        fi

        # if crams_with_multiple_samples is not empty, write an error message to qc_messages.txt
        n_multi_sample_crams=${#crams_with_multiple_samples[@]}
        if [ $n_multi_sample_crams -ne 0 ]; then
            {
                if [ $n_multi_sample_crams -eq 1 ]; then
                    pluralized=""
                else
                    pluralized="s"
                fi
                joined=$(IFS=","; echo "${crams_with_multiple_samples[*]}")
                list_to_show="${joined//,/, }" # Replaces every ',' with ', '
                echo "Found $n_multi_sample_crams CRAM file$pluralized containing data for more than one sample: $list_to_show"
            } >> qc_messages.txt
        else
            echo "All CRAM files contain data for exactly one sample."
        fi

        # if duplicate_sample_id_descriptions is not empty, write an error message to qc_messages.txt
        n_duplicate_sample_ids=${#duplicate_sample_id_descriptions[@]}
        if [ $n_duplicate_sample_ids -ne 0 ]; then
            {
                if [ $n_duplicate_sample_ids -eq 1 ]; then
                    pluralized=""
                else
                    pluralized="s"
                fi
                # each description already contains its own ", "-separated CRAM list, so join with a
                # literal ", " directly rather than the "replace every comma" trick used above (that
                # would also mangle the commas inside each description's CRAM list)
                list_to_show=$(printf '%s, ' "${duplicate_sample_id_descriptions[@]}")
                list_to_show="${list_to_show%, }" # strip the trailing ", " left by the printf loop
                echo "Found $n_duplicate_sample_ids duplicate sample ID$pluralized across CRAMs: $list_to_show"
            } >> qc_messages.txt
        else
            echo "All CRAM sample IDs are unique across the provided CRAMs."
        fi

        # passes_qc is true if qc_messages is empty
        if [ ! -s qc_messages.txt ]; then
            echo "true" > passes_qc.txt
        else
            echo "false" > passes_qc.txt
        fi

        # This task should always succeed
        exit 0
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23.1"
        cpu: cpu
        disks: "local-disk 10 HDD"
        memory: "4 GiB"
        maxRetries: 2
    }

    output {
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}
