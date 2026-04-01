version 1.0

workflow InputQC {
    input {
        String pipeline_version = "0.0.1"

        Array[String] contigs

        # this is the path the a directory that contains sites vcf, sites tabke, and reference chunks file.  should end with a "/
        String reference_panel_prefix

        # service currently does not accept VCFs as input
        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? cram_manifest
        File fasta
        File fasta_index
        String output_basename

        File ref_dict
    }

    # validate that either crams, or cram manifest is provided
    if (!defined(crams) && !defined(cram_manifest)) {
        Boolean no_data_passes_qc = false
        String no_data_message = "No input data provided. Please provide either CRAM files or a CRAM manifest."
    }

    # validate that not more than one of these is provided
    if ( (defined(crams) && defined(cram_manifest))) {
        Boolean multiple_data_types_passes_qc = false
        String multiple_data_types_message = "Multiple input data types provided. Please provide only CRAM files (with corresponding CRAIs and sample IDs) or a CRAM manifest."
    }

    # validations for array crams input
    if (defined(crams) && !defined(cram_manifest)) {
        if (!defined(cram_indices) || !defined(sample_ids)) {
            Boolean no_cram_index_or_sample_id_passes_qc = false
            String no_cram_index_or_sample_id_message = "CRAM indices and sample IDs are required when CRAM files are provided. Please provide cram index files and a list of sample IDs corresponding to the CRAM files."
        }

        if (defined(crams) && defined(cram_indices) && defined(sample_ids)) {
            call ValidateCramsAndIndices {
                input:
                    crams = select_first([crams]),
                    cram_indices = select_first([cram_indices]),
                    sample_ids = select_first([sample_ids])
            }
        }
        
    }

    # validations for cram manifest input
    if (defined(cram_manifest)) {
        call ValidateCramManifest {
            input:
                cram_manifest = select_first([cram_manifest])
        }
    }

    output {
        Boolean passes_qc = select_first([no_data_passes_qc, multiple_data_types_passes_qc, no_cram_index_or_sample_id_passes_qc, ValidateCramsAndIndices.passes_qc, ValidateCramManifest.passes_qc])
        String qc_messages = select_first([no_data_message, multiple_data_types_message, no_cram_index_or_sample_id_message, ValidateCramsAndIndices.qc_messages, ValidateCramManifest.qc_messages])
    }
}


task ValidateCramsAndIndices {
    input {
        Array[File] crams
        Array[File] cram_indices
        Array[String] sample_ids

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(1.1*size(crams, "GiB")) + 10
    }

    Int num_crams = length(crams)
    Int num_cram_indices = length(cram_indices)
    Int num_sample_ids = length(sample_ids)

    command <<<
        # create empty qc messages file
        touch qc_messages.txt

        # validate that the number of CRAMs, CRAIs, and sample IDs match
        if [ ~{num_crams} -ne ~{num_cram_indices} ] || [ ~{num_crams} -ne ~{num_sample_ids} ]; then
            echo "Mismatch in the number of CRAMs (~{num_crams}), CRAIs (~{num_cram_indices}), and sample IDs (~{num_sample_ids})." >> qc_messages.txt
        else
            echo "Number of CRAMs, CRAIs, and sample IDs match: ~{num_crams}."
        fi

        # validate that sample IDs are unique
        unique_sample_ids=$(printf "%s\n" "${sample_ids[@]}" | sort -u | wc -l)
        if [ $unique_sample_ids -ne ~{num_sample_ids} ]; then
            # find duplicate sample IDs
            duplicate_sample_ids=$(printf "%s\n" "${sample_ids[@]}" | sort | uniq -d)
            echo "Duplicate sample IDs found: ${duplicate_sample_ids}" >> qc_messages.txt
        else
            echo "Sample IDs are unique."
        fi

        # ensure all crams end with .cram and all cram indices end with .crai
        crams_with_wrong_extension=$(printf "%s\n" "${crams[@]}" | grep -vE "\.cram$")
        if [ ! -z "${crams_with_wrong_extension}" ]; then
            echo "The following CRAM files do not have a .cram extension: ${crams_with_wrong_extension}" >> qc_messages.txt
        else
            echo "All CRAM files have the correct .cram extension."
        fi
        cram_indices_with_wrong_extension=$(printf "%s\n" "${cram_indices[@]}" | grep -vE "\.crai$")
        if [ ! -z "${cram_indices_with_wrong_extension}" ]; then
            echo "The following CRAM index files do not have a .crai extension: ${cram_indices_with_wrong_extension}" >> qc_messages.txt
        else
            echo "All CRAM index files have the correct .crai extension."
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
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 3
        maxRetries: 1
        noAddress: true
    }
    
    output {
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}

task ValidateCramManifest {
    input {
        File cram_manifest

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(1.1*size(cram_manifest, "GiB")) + 10
    }

    command <<<
        # create empty qc messages file
        touch qc_messages.txt

        # validate that the manifest has the required columns, independent of order
        required_columns=("sample_id" "cram_path" "cram_index_path")
        header=$(head -n 1 ${cram_manifest})
        missing_columns=()

        for col in "${required_columns[@]}"; do
            if ! echo "${header}" | tr '\t' '\n' | grep -q "^${col}$"; then
                missing_columns+=("${col}")
            fi
        done

        if [ ${#missing_columns[@]} -gt 0 ]; then
            echo "CRAM manifest is missing required columns: ${missing_columns[*]}. Expected columns: ${required_columns[*]}. Actual header: ${header}" >> qc_messages.txt
        else
            echo "CRAM manifest has the required columns."
        fi

        # validate that sample IDs are unique
        # find the column number for sample_id
        sample_id_col=$(head -n 1 ${cram_manifest} | tr '\t' '\n' | grep -n "^sample_id$" | cut -d: -f1)
        
        if [ -z "${sample_id_col}" ]; then
            echo "Unable to determine sample_id column position." >> qc_messages.txt
        else
            sample_ids=$(tail -n +2 ${cram_manifest} | cut -f${sample_id_col})
            unique_sample_ids=$(echo "${sample_ids}" | sort -u | wc -l)
            total_sample_ids=$(echo "${sample_ids}" | wc -l)
            if [ $unique_sample_ids -ne $total_sample_ids ]; then
                # find duplicate sample IDs
                duplicate_sample_ids=$(echo "${sample_ids}" | sort | uniq -d)
                echo "Duplicate sample IDs found in CRAM manifest: ${duplicate_sample_ids}" >> qc_messages.txt
            else
                echo "Sample IDs in CRAM manifest are unique."
            fi
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
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 3
        maxRetries: 1
        noAddress: true
    }
    
    output {
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}
