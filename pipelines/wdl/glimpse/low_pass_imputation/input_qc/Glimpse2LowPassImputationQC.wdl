version 1.0

workflow InputQC {
    # if this changes, update the input_qc_version value in Glimpse2LowPassImputation.wdl
    String pipeline_version = "1.0.2"

    input {
        # service expects only cram_manifest even though main wdl can alternatively take input arrays
        File cram_manifest
        String output_basename

        Array[String] contigs
        # this is the path to a directory that contains sites vcf, sites table, and reference chunks file. should end with a "/"
        String reference_panel_prefix
        File fasta
        File fasta_index
        File ref_dict

        # used for warp tests only (which use inputs in an RP bucket). service does not support RP buckets and will not provide this input.
        String? billing_project_for_rp
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
                sample_ids = ConvertCramManifestToInputArrays.sample_ids,
                billing_project_for_rp = billing_project_for_rp
        }
    }

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
        sample_ids_filename = "sample_ids.txt"

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
            required_cols = ['sample_id', 'cram_path', 'cram_index_path']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                with open(qc_messages_filename, 'w') as qc_file:
                    qc_file.write(f"Missing required columns in the CRAM manifest: {', '.join(missing_cols)}.")
                with open(passes_qc_filename, 'w') as f:
                    f.write("false")
                
                # Create empty output files
                open(crams_filename, 'w').close()
                open(cram_indices_filename, 'w').close()
                open(sample_ids_filename, 'w').close()
            else:
                # Write to output files, stripping leading/trailing whitespace from each value
                write_column(df['sample_id'], sample_ids_filename)
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
            open(sample_ids_filename, 'w').close()

        EOF
        python3 script.py

        # This task should always succeed
        exit 0
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
        Array[String] sample_ids = read_lines("sample_ids.txt")
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}


task ValidateCramsAndIndicesAndSampleIds {
    input {
        Array[String] crams
        Array[String] cram_indices
        Array[String] sample_ids

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
        parsed_sample_ids = """~{sep=' ' sample_ids}""".split()
        parsed_crams = """~{sep=' ' crams}""".split()
        parsed_cram_indices = """~{sep=' ' cram_indices}""".split()

        # remove empty strings
        sample_ids = [s for s in parsed_sample_ids if s]
        crams = [c for c in parsed_crams if c]
        cram_indices = [c for c in parsed_cram_indices if c]

        num_crams = len(crams)
        num_cram_indices = len(cram_indices)
        num_sample_ids = len(sample_ids)

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

        # Validate that the number of CRAMs, CRAIs, and sample IDs match
        if num_crams != num_cram_indices or num_crams != num_sample_ids:
            qc_messages.append(f"Found different numbers of CRAMs ({num_crams}), CRAM index files ({num_cram_indices}), and sample IDs ({num_sample_ids}).")
        else:
            print(f"Number of CRAMs, CRAM index files, and sample IDs match: found {num_crams} of each.")

        # Validate that sample IDs are unique
        unique_sample_ids = set(sample_ids)
        if len(unique_sample_ids) != num_sample_ids:
            duplicates = [sid for sid in unique_sample_ids if sample_ids.count(sid) > 1]
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(duplicates), 'duplicate sample ID')}", duplicates))
        else:
            print("Sample IDs are unique.")

        # Ensure all crams end with .cram and all cram indices end with .crai
        crams_with_wrong_extension = [c for c in crams if not c.endswith('.cram')]
        if crams_with_wrong_extension:
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(crams_with_wrong_extension), 'CRAM file')} that do not have a .cram extension", crams_with_wrong_extension))
        else:
            print("All CRAM files have the correct .cram extension.")

        cram_indices_with_wrong_extension = [c for c in cram_indices if not c.endswith('.crai')]
        if cram_indices_with_wrong_extension:
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(cram_indices_with_wrong_extension), 'CRAM index file')} without a .crai extension", cram_indices_with_wrong_extension))
        else:
            print("All CRAM index files have the correct .crai extension.")

        # Validate that cram paths are unique
        unique_crams = set(crams)
        if len(unique_crams) != num_crams:
            duplicates = [c for c in unique_crams if crams.count(c) > 1]
            qc_messages.append(create_error_message_with_item_list(f"Found {pluralize(len(duplicates), 'set')} of duplicate CRAM paths", duplicates))
        else:
            print("CRAM paths are unique.")

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

        # This task should always succeed
        exit 0
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
    }

    String billing_project = select_first([billing_project_for_rp, ""])

    command <<<
        # set up auth for accessing files using samtools
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # configure billing project to use for requester pays buckets, if billing project provided
        if [ -n "$billing_project" ]; then
            echo "Using billing project '$billing_project' for requester pays buckets."
            export GCS_REQUESTER_PAYS_PROJECT=~{billing_project}
        fi

        # collect reference MD5sums from the reference dictionary for the contigs in the reference panel
        contigs=('chr1' 'chr2')
        ref_dict='/scripts/Homo_sapiens_assembly38.dict'
        cram='/scripts/NA12878_PLUMBING.cram'


        declare -A ref_md5sums
        while read -r line; do
            if [[ $line == @SQ* ]]; then
                chrom=$(echo "$line" | awk '{for(i=1;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}')
                md5=$(echo "$line" | awk '{for(i=1;i<=NF;i++) if($i~/^M5:/) print substr($i,4)}')

                if [[ " ${contigs[@]} " =~ " ${chrom} " ]]; then
                    ref_md5sums["$chrom"]="$md5"

                    echo  ${contigs[@]}
                    echo ${chrom}
                    
                    # Check if we've found all contigs
                    if [[ ${#ref_md5sums[@]} -eq ${#contigs[@]} ]]; then
                        break
                    fi
                fi
            fi
        done < $ref_dict

        echo "found relevant contigs with these md5sums in ref dict ${ref_dict}:"
        for chrom in "${!ref_md5sums[@]}"; do
            echo "  $chrom: ${ref_md5sums[$chrom]}"
        done

        # read cram headers to validate that they contain the expected reference alignment MD5sums
        echo "Validating CRAM file: $cram"
        header=$(samtools view -H "$cram")
        for chrom in "${!ref_md5sums[@]}"; do
            expected_md5=${ref_md5sums[$chrom]}
            if ! echo "$header" | grep -q "SN:$chrom.*M5:$expected_md5"; then
                echo "CRAM file $cram is missing expected reference alignment MD5 for contig $chrom (expected SN:$chrom and M5:$expected_md5 on the same line in header)" >> qc_messages.txt
            else
                echo "CRAM file $cram contains expected reference alignment MD5 for contig $chrom."
            fi
        done

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
        docker: "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
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
