version 1.0

workflow InputQC {
    # if this changes, update the input_qc_version value in Glimpse2LowPassImputation.wdl
    String pipeline_version = "0.0.1"

    input {
        Array[String] contigs

        # this is the directory that contains sites vcf, sites table, and reference chunks file. should end with a "/"
        String reference_panel_prefix

        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? cram_manifest
        String output_basename

        File fasta
        File fasta_index
        File ref_dict

        String? billing_project_for_rp
    }

    # validate that either crams, or cram manifest is provided
    if (!defined(crams) && !defined(cram_manifest)) {
        Boolean no_data_passes_qc = false
        String no_data_message = "No input data provided. Please provide either CRAM files or a CRAM manifest."
    }

    # validate that not more than one of these is provided
    if ((defined(crams) && defined(cram_manifest))) {
        Boolean multiple_data_types_passes_qc = false
        String multiple_data_types_message = "Multiple input data types provided. Please provide only CRAM files (with corresponding CRAM index files and sample IDs) or a CRAM manifest."
    }

    # convert cram manifest to arrays of crams, cram indices, and sample ids if manifest is provided
    if (defined(cram_manifest) && !defined(crams)) {
        call ConvertCramManifestToInputArrays {
            input:
                cram_manifest = select_first([cram_manifest])
        }
    }

    Array[String] cram_array = select_first([crams, ConvertCramManifestToInputArrays.crams, []])
    Array[String] cram_indices_array = select_first([cram_indices, ConvertCramManifestToInputArrays.cram_indices, []])
    Array[String] sample_ids_array = select_first([sample_ids, ConvertCramManifestToInputArrays.sample_ids, []])

    if (length(cram_array) > 0) {
        call ValidateCramsAndIndicesAndSampleIds {
            input:
                crams = cram_array,
                cram_indices = cram_indices_array,
                sample_ids = sample_ids_array,
                billing_project_for_rp = billing_project_for_rp
        }
        # Boolean do_cram_qc = select_first([ConvertCramManifestToInputArrays.passes_qc, defined(crams) && !defined(cram_manifest), false])
        
        # validations for array crams, cram indices, and sample ids (whether supplied directly or via manifest)
        # if (do_cram_qc) {
            # Array[String] cram_array = select_first([crams, ConvertCramManifestToInputArrays.crams])

            # Array[String] evaluated_cram_indices = select_first([cram_indices, ConvertCramManifestToInputArrays.cram_indices, []])
            # Array[String] evaluated_sample_ids = select_first([sample_ids, ConvertCramManifestToInputArrays.sample_ids, []])

            # Boolean cram_indices_and_sample_ids_provided = (length(evaluated_cram_indices) > 0) && (length(evaluated_sample_ids) > 0)
            
            # if (!cram_indices_and_sample_ids_provided) {
            #     Boolean no_cram_index_or_sample_id_passes_qc = false
            #     String no_cram_index_or_sample_id_message = "CRAM indices and sample IDs are required when CRAM files are provided. Please provide both CRAM index files and a corresponding list of sample IDs."
            # }

            # if (cram_indices_and_sample_ids_provided) {
                # call ValidateCramsAndIndicesAndSampleIds {
                #     input:
                #         crams = cram_array,
                #         cram_indices = evaluated_cram_indices,
                #         sample_ids = evaluated_sample_ids,
                #         billing_project_for_rp = billing_project_for_rp
                # }
            # }
        # }
    }

    output {
        Boolean passes_qc = select_first([ValidateCramsAndIndicesAndSampleIds.passes_qc, ConvertCramManifestToInputArrays.passes_qc, no_data_passes_qc, multiple_data_types_passes_qc])
        String qc_messages = select_first([ValidateCramsAndIndicesAndSampleIds.qc_messages, ConvertCramManifestToInputArrays.qc_messages, no_data_message, multiple_data_types_message])
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

        # Read the manifest
        try:
            df = pd.read_csv("~{cram_manifest}", sep='\t')
            
            # Check for required columns
            required_cols = ['sample_id', 'cram_path', 'cram_index_path']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                with open(qc_messages_filename, 'w') as qc_file:
                    qc_file.write(f"Missing required columns in the CRAM manifest: {', '.join(missing_cols)}")
                with open(passes_qc_filename, 'w') as f:
                    f.write("false")
                
                # Create empty output files
                open(crams_filename, 'w').close()
                open(cram_indices_filename, 'w').close()
                open(sample_ids_filename, 'w').close()
            else:
                # Write to output files
                df['sample_id'].to_csv(sample_ids_filename, index=False, header=False)
                df['cram_path'].to_csv(crams_filename, index=False, header=False)
                df['cram_index_path'].to_csv(cram_indices_filename, index=False, header=False)
                
                # Write QC results
                with open(qc_messages_filename, 'w') as f:
                    f.write('\n'.join(qc_messages) if qc_messages else '')
                
                with open(passes_qc_filename, 'w') as f:
                    f.write("true" if not qc_messages else "false")
        
        except Exception as e:
            with open(qc_messages_filename, 'w') as qc_file:
                qc_file.write(f"Error reading CRAM manifest: {str(e)}")
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
        memory: "1 GiB"
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

    command <<<
        pip install google-cloud-storage

        cat <<'EOF' > script.py
        from google.cloud import storage

        qc_messages = []

        # Parse WDL arrays from space-separated strings
        sample_ids = """~{sep=' ' sample_ids}""".split()
        crams = """~{sep=' ' crams}""".split()
        cram_indices = """~{sep=' ' cram_indices}""".split()

        num_crams = len(crams)
        num_cram_indices = len(cram_indices)
        num_sample_ids = len(sample_ids)

        # Validate that the number of CRAMs, CRAIs, and sample IDs match
        if num_crams != num_cram_indices or num_crams != num_sample_ids:
            qc_messages.append(f"Found different numbers of CRAMs ({num_crams}), CRAIs ({num_cram_indices}), and sample IDs ({num_sample_ids}).")
        else:
            print(f"Number of CRAMs, CRAIs, and sample IDs match: found {num_crams} of each.")

        # Validate that sample IDs are unique
        unique_sample_ids = set(sample_ids)
        if len(unique_sample_ids) != num_sample_ids:
            duplicates = [sid for sid in unique_sample_ids if sample_ids.count(sid) > 1]
            qc_messages.append(f"Duplicate sample IDs found: {', '.join(sorted(duplicates))}")
        else:
            print("Sample IDs are unique.")

        # Ensure all crams end with .cram and all cram indices end with .crai
        crams_with_wrong_extension = [c for c in crams if not c.endswith('.cram')]
        if crams_with_wrong_extension:
            qc_messages.append(f"The following CRAM files do not have a .cram extension: {', '.join(crams_with_wrong_extension)}")
        else:
            print("All CRAM files have the correct .cram extension.")

        cram_indices_with_wrong_extension = [c for c in cram_indices if not c.endswith('.crai')]
        if cram_indices_with_wrong_extension:
            qc_messages.append(f"The following CRAM index files do not have a .crai extension: {', '.join(cram_indices_with_wrong_extension)}")
        else:
            print("All CRAM index files have the correct .crai extension.")

        # Validate that cram paths are unique
        unique_crams = set(crams)
        if len(unique_crams) != num_crams:
            duplicates = [c for c in unique_crams if crams.count(c) > 1]
            qc_messages.append(f"Duplicate CRAM paths found: {', '.join(sorted(duplicates))}")
        else:
            print("CRAM paths are unique.")

        # Ensure that all CRAM files are less than the maximum file size allowed
        max_cram_file_size_gb = ~{max_cram_file_size_gb}
        billing_project = "~{billing_project_for_rp}"
        crams_exceeding_max_size = []

        client = storage.Client(project=billing_project) if billing_project else storage.Client()
        for cram in crams:
            try:
                if cram.startswith('gs://'):
                    blob = storage.Blob.from_uri(cram, client=client)
                    
                    # Reload to get metadata
                    # blob.reload()
                    
                    # Get file size
                    file_size_bytes = blob.size
                    file_size_gb = file_size_bytes // (1024 ** 3)
                    print(f" - File size for {cram}: {file_size_gb} GB")
                    
                    if file_size_gb > max_cram_file_size_gb:
                        crams_exceeding_max_size.append(f"{cram} ({file_size_gb}GB)")
                else:
                    qc_messages.append(f"Invalid GCS path format for {cram}. Expected gs:// prefix.")
            except Exception as e:
                qc_messages.append(f"Error checking file size for {cram}: {str(e)}")

        if crams_exceeding_max_size:
            qc_messages.append(f"The following CRAM files exceed the maximum allowed file size of {max_cram_file_size_gb}GB: {', '.join(crams_exceeding_max_size)}")
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
        memory: "1 GiB"
        preemptible: 3
        maxRetries: 2
    }
    
    output {
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}
