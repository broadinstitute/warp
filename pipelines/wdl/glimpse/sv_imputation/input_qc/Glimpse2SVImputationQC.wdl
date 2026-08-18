version 1.0

workflow InputQC {
    # if this changes, update the input_qc_version value in Glimpse2SVImputation.wdl
    String pipeline_version = "0.0.1"

    input {
        # service expects only gvcf_manifest even though main wdl can alternatively take input arrays
        File gvcf_manifest
        String output_prefix

        # remaining inputs kept for interface consistency with Glimpse2SVImputation.wdl; not all are used by QC
        File preprocess_panel_bubble_split_sites_only_vcf
        File preprocess_panel_bubble_split_sites_only_vcf_idx

        Array[String] paste_regions

        Array[String] chromosomes
        File genetic_maps_tsv
        File ref_dict
        File chunked_panel_json

        File pop_glimpse2_panel_resources_json

        Float? info_filter_for_inclusion

        # used for warp tests only (which use inputs in an RP bucket). service does not support RP buckets and will not provide this input.
        String? billing_project_for_rp
        Int? max_file_size_gb_override

        # optional additional header line to add to the output VCF
        String? pipeline_header_line
    }

    call ValidateGvcfManifest {
        input:
            gvcf_manifest = gvcf_manifest,
            billing_project_for_rp = billing_project_for_rp
            max_gvcf_file_size_gb = max_file_size_gb_override
    }

    # only validate individual GVCF contents if the manifest itself passed QC
    if (ValidateGvcfManifest.passes_qc) {
        call ValidateGvcfInput {
            input:
                gvcfs = ValidateGvcfManifest.gvcfs,
                ref_dict = ref_dict,
                billing_project_for_rp = billing_project_for_rp
        }
    }

    output {
        Boolean passes_qc = select_first([ValidateGvcfInput.passes_qc, ValidateGvcfManifest.passes_qc])
        String qc_messages = select_first([ValidateGvcfInput.qc_messages, ValidateGvcfManifest.qc_messages])
    }
}


task ValidateGvcfManifest {
    input {
        File gvcf_manifest

        Int max_gvcf_file_size_gb = 50
        String? billing_project_for_rp # if set, will use this to check file sizes for requester pays buckets. if not set and input is in a RP bucket, the check will fail
    }

    String billing_project = select_first([billing_project_for_rp, ""])

    command <<<
        pip install google-cloud-storage

        cat <<EOF > script.py
        import pandas as pd
        from google.cloud import storage

        max_gvcf_file_size_gb = ~{max_gvcf_file_size_gb}
        billing_project = "~{billing_project}"

        qc_messages = []

        qc_messages_filename = "qc_messages.txt"
        passes_qc_filename = "passes_qc.txt"
        gvcfs_filename = "gvcfs.txt"
        gvcf_indices_filename = "gvcf_indices.txt"

        MAX_ITEMS_IN_ERROR_MESSAGES = 5

        def create_error_message_with_item_list(base_error_message, items_list):
            """Helper function to create error messages that include a list of items, but truncates the list if it's too long."""
            items_list_to_show = items_list[:MAX_ITEMS_IN_ERROR_MESSAGES]
            exceeded_limit_message = f"; first {MAX_ITEMS_IN_ERROR_MESSAGES} are" if len(items_list) > MAX_ITEMS_IN_ERROR_MESSAGES else ""
            return f"{base_error_message}{exceeded_limit_message}: {', '.join(items_list_to_show)}."

        def pluralize(number, subject):
            """Helper function to return a properly pluralized phrase based on the number provided, e.g. '1 GVCF file' or '2 GVCF files'."""
            return f"{number} {subject}" if number == 1 else f"{number} {subject}s"

        def write_column(column_data, filename):
            """Write column to file, with each value stripped of leading/trailing whitespace."""
            filtered = column_data.fillna('').astype(str).str.strip()
            with open(filename, 'w') as f:
                for value in filtered:
                    f.write(f"{value}\n")

        def write_outputs(qc_messages, gvcfs=None, gvcf_indices=None):
            with open(qc_messages_filename, 'w') as f:
                f.write('\n'.join(qc_messages) if qc_messages else '')
            with open(passes_qc_filename, 'w') as f:
                f.write("true" if not qc_messages else "false")
            if gvcfs is not None:
                write_column(gvcfs, gvcfs_filename)
            else:
                open(gvcfs_filename, 'w').close()
            if gvcf_indices is not None:
                write_column(gvcf_indices, gvcf_indices_filename)
            else:
                open(gvcf_indices_filename, 'w').close()

        try:
            df = pd.read_csv("~{gvcf_manifest}", sep='\t')

            required_cols = ['gvcf_path', 'gvcf_index_path']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                write_outputs([f"Missing required column header(s) in the GVCF manifest: {', '.join(missing_cols)}."])
            elif len(df) == 0:
                write_outputs(["The GVCF manifest must contain at least one row."])
            elif df[required_cols].isnull().any().any():
                write_outputs(["The GVCF manifest contains empty values in required columns."])
            else:
                gvcf_paths = df['gvcf_path'].astype(str).str.strip()
                gvcf_index_paths = df['gvcf_index_path'].astype(str).str.strip()

                # Validate that GVCF paths are unique
                duplicate_gvcf_paths = [g for g in set(gvcf_paths) if list(gvcf_paths).count(g) > 1]
                if duplicate_gvcf_paths:
                    qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(duplicate_gvcf_paths), 'set')} of duplicate GVCF paths", duplicate_gvcf_paths))
                else:
                    print("GVCF paths are unique.")

                # Ensure all GVCFs and indices have the expected extensions
                gvcfs_with_wrong_extension = [g for g in gvcf_paths if not g.endswith('.vcf.gz')]
                if gvcfs_with_wrong_extension:
                    qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(gvcfs_with_wrong_extension), 'GVCF file')} without a .vcf.gz extension", gvcfs_with_wrong_extension))
                else:
                    print("All GVCF files have the correct .vcf.gz extension.")

                gvcf_indices_with_wrong_extension = [g for g in gvcf_index_paths if not g.endswith('.tbi')]
                if gvcf_indices_with_wrong_extension:
                    qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(gvcf_indices_with_wrong_extension), 'GVCF index file')} without a .tbi extension", gvcf_indices_with_wrong_extension))
                else:
                    print("All GVCF index files have the correct .tbi extension.")

                # Validate that each gvcf-index pair has matching basenames
                mismatched_basename_pairs = []
                for gvcf, index in zip(gvcf_paths, gvcf_index_paths):
                    gvcf_basename = gvcf.split('/')[-1]
                    index_basename = index.split('/')[-1]
                    if not index_basename.startswith(gvcf_basename):
                        mismatched_basename_pairs.append(f"{gvcf} and {index}")
                if mismatched_basename_pairs:
                    qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(mismatched_basename_pairs), 'GVCF-index pair')} with mismatched basenames", mismatched_basename_pairs))
                else:
                    print("All GVCF-index pairs have matching basenames.")

                # Ensure that all GVCF and index files exist and are accessible, and that GVCFs are within the max file size
                files_with_invalid_gcs_format = []
                files_with_access_issues = []
                gvcfs_exceeding_max_size = []

                client = storage.Client()

                def check_blob(path):
                    if not path.startswith('gs://'):
                        files_with_invalid_gcs_format.append(path)
                        return None
                    bucket_name, blob_name = path[5:].split('/', 1)
                    bucket = client.bucket(bucket_name, user_project=billing_project) if billing_project else client.bucket(bucket_name)
                    blob = bucket.blob(blob_name)
                    try:
                        blob.reload(client=client)
                    except Exception as e:
                        files_with_access_issues.append(path)
                        print(f"ERROR DETAILS for {path}: {str(e)}")
                        return None
                    return blob

                for gvcf_path, gvcf_index_path in zip(gvcf_paths, gvcf_index_paths):
                    gvcf_blob = check_blob(gvcf_path)
                    check_blob(gvcf_index_path)

                    if gvcf_blob:
                        file_size_gb = int(gvcf_blob.size) / (1024 ** 3)
                        print(f"File size for {gvcf_path}: {file_size_gb:.2f} GB")
                        if file_size_gb > max_gvcf_file_size_gb:
                            gvcfs_exceeding_max_size.append(f"{gvcf_path} ({file_size_gb:.2f}GB)")

                if files_with_invalid_gcs_format:
                    qc_messages.append(create_error_message_with_item_list(
                        f"Found {pluralize(len(files_with_invalid_gcs_format), 'file')} with invalid GCS format (must start with 'gs://')",
                        files_with_invalid_gcs_format))

                if files_with_access_issues:
                    qc_messages.append(create_error_message_with_item_list(
                        f"Found {pluralize(len(files_with_access_issues), 'file')} that could not be accessed (may be due to non-existent files, lack of permissions, or requester pays bucket)",
                        files_with_access_issues))

                if gvcfs_exceeding_max_size:
                    qc_messages.append(create_error_message_with_item_list(
                        f"Found {pluralize(len(gvcfs_exceeding_max_size), 'GVCF file')} exceeding the maximum allowed file size of {max_gvcf_file_size_gb}GB",
                        gvcfs_exceeding_max_size))
                else:
                    print(f"All GVCF files are within the maximum allowed file size of {max_gvcf_file_size_gb}GB.")

                write_outputs(qc_messages, gvcf_paths, gvcf_index_paths)

        except Exception as e:
            write_outputs([f"Error reading GVCF manifest: {str(e)}."])

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
    }

    output {
        Array[String] gvcfs = read_lines("gvcfs.txt")
        Array[String] gvcf_indices = read_lines("gvcf_indices.txt")
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}


task ValidateGvcfInput {
    input {
        Array[String] gvcfs
        File ref_dict

        String? billing_project_for_rp # if set, will use this to access GVCFs in requester pays buckets. if not set and input is in a RP bucket, the check will fail
    }

    String billing_project = select_first([billing_project_for_rp, ""])
    String ref_dict_basename = basename(ref_dict)

    command <<<
        set -uo pipefail

        # set up auth for accessing files using bcftools
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # configure billing project to use for requester pays buckets, if billing project provided
        if [ -n "~{billing_project}" ]; then
            echo "Using billing project '~{billing_project}' for requester pays buckets."
            export GCS_REQUESTER_PAYS_PROJECT=~{billing_project}
        fi

        touch qc_messages.txt

        ref_dict_basename="~{ref_dict_basename}"

        MAX_ITEMS_IN_ERROR_MESSAGES=5
        # to limit runtime of this task, we will only check the first MAX_GVCFS_TO_CHECK GVCFs
        MAX_GVCFS_TO_CHECK=100

        # Appends a truncated, comma-separated summary of $2... to qc_messages.txt, prefixed by $1, if any items are given.
        append_aggregated_message() {
            local base_message="$1"
            shift
            local items=("$@")
            local n_items=${#items[@]}
            if [ "$n_items" -eq 0 ]; then
                return
            fi
            local joined
            if [ "$n_items" -gt "$MAX_ITEMS_IN_ERROR_MESSAGES" ]; then
                joined=$(IFS=","; echo "${items[*]:0:$MAX_ITEMS_IN_ERROR_MESSAGES}")
                echo "${base_message}; first ${MAX_ITEMS_IN_ERROR_MESSAGES} are: ${joined//,/, }." >> qc_messages.txt
            else
                joined=$(IFS=","; echo "${items[*]}")
                echo "${base_message}: ${joined//,/, }." >> qc_messages.txt
            fi
        }

        gvcfs_with_incompatible_contigs=()
        gvcfs_with_missing_format_fields=()
        gvcf_check_count=0

        # this task assumes the manifest check has already confirmed that all GVCFs exist and are accessible
        for gvcf in ~{sep=' ' gvcfs}; do
            gvcf_check_count=$((gvcf_check_count + 1))
            echo "Validating GVCF file: $gvcf"

            # --validation-type-to-exclude ALL skips variant-level validation and only checks that the
            # VCF header's sequence dictionary is compatible with the provided reference dictionary
            gatk ValidateVariants \
                -V "$gvcf" \
                --validate-GVCF true \
                --sequence-dictionary ~{ref_dict} \
                --validation-type-to-exclude ALL \
                ~{"--gcs-project-for-requester-pays " + billing_project_for_rp} \
                2> gatk_output.txt

            cat gatk_output.txt

            if grep -q "incompatible contigs" gatk_output.txt; then
                echo "GVCF file $gvcf has contigs incompatible with the expected reference dictionary ($ref_dict_basename)."
                gvcfs_with_incompatible_contigs+=("$gvcf")
            else
                echo "GVCF file $gvcf has contigs compatible with the expected reference dictionary."
            fi

            # Ensure the PL and GT FORMAT annotations are declared in the header.
            header=$(bcftools view -h "$gvcf")
            format_lines=$(echo "$header" | grep '^##FORMAT=<')
            missing_format_fields=()
            if ! echo "$format_lines" | grep -q 'ID=PL[,>]'; then
                missing_format_fields+=("PL")
            fi
            if ! echo "$format_lines" | grep -q 'ID=GT[,>]'; then
                missing_format_fields+=("GT")
            fi
            if [ ${#missing_format_fields[@]} -gt 0 ]; then
                echo "GVCF file $gvcf is missing expected FORMAT annotation(s) in its header: ${missing_format_fields[*]}"
                gvcfs_with_missing_format_fields+=("$gvcf")
            else
                echo "GVCF file $gvcf declares the expected PL and GT FORMAT annotations in its header."
            fi

            # if we've found more than MAX_ITEMS_IN_ERROR_MESSAGES gvcfs with issues, we can stop checking the rest of the gvcfs because the error messages will be truncated anyway
            total_issue_count=$(( ${#gvcfs_with_incompatible_contigs[@]} + ${#gvcfs_with_missing_format_fields[@]} ))
            if [ "$total_issue_count" -gt "$MAX_ITEMS_IN_ERROR_MESSAGES" ]; then
                echo "Found more than $MAX_ITEMS_IN_ERROR_MESSAGES GVCF files with issues; skipping validation of remaining GVCF files"
                break
            fi
            # if we've checked more than MAX_GVCFS_TO_CHECK gvcfs, we will stop to limit runtime of this task
            if [ "$gvcf_check_count" -ge "$MAX_GVCFS_TO_CHECK" ]; then
                echo "Checked $MAX_GVCFS_TO_CHECK GVCF files; stopping further checks to limit runtime of this task"
                break
            fi
        done

        n_bad_contig_gvcfs=${#gvcfs_with_incompatible_contigs[@]}
        if [ "$n_bad_contig_gvcfs" -ne 0 ]; then
            if [ "$n_bad_contig_gvcfs" -eq 1 ]; then pluralized=""; else pluralized="s"; fi
            append_aggregated_message "Found $n_bad_contig_gvcfs GVCF file$pluralized with contigs incompatible with the expected reference dictionary ($ref_dict_basename)" "${gvcfs_with_incompatible_contigs[@]}"
        else
            echo "All checked GVCF files have contigs compatible with the expected reference dictionary."
        fi

        n_missing_format_gvcfs=${#gvcfs_with_missing_format_fields[@]}
        if [ "$n_missing_format_gvcfs" -ne 0 ]; then
            if [ "$n_missing_format_gvcfs" -eq 1 ]; then pluralized=""; else pluralized="s"; fi
            append_aggregated_message "Found $n_missing_format_gvcfs GVCF file$pluralized missing the required PL and/or GT FORMAT annotation(s) in its header" "${gvcfs_with_missing_format_fields[@]}"
        else
            echo "All checked GVCF files declare the expected PL and GT FORMAT annotations in their headers."
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
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "4 GiB"
        maxRetries: 2
        noAddress: true
    }

    output {
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}
