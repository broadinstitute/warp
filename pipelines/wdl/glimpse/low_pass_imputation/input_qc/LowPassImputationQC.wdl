version 1.0

workflow InputQC {
    # if this changes, update the input_qc_version value in Glimpse2LowPassImputation.wdl
    String pipeline_version = "0.0.1"

    input {

        Array[String] contigs

        # this is the path the a directory that contains sites vcf, sites tabke, and reference chunks file.  should end with a "/"
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
    if ((defined(crams) && defined(cram_manifest))) {
        Boolean multiple_data_types_passes_qc = false
        String multiple_data_types_message = "Multiple input data types provided. Please provide only CRAM files (with corresponding CRAM index files and sample IDs) or a CRAM manifest."
    }

    # convert cram manifest to arrays of crams, cram indices, and sample ids if manifest is provided
    if (defined(cram_manifest)) {
        call ConvertCramManifestToCramArrays {
            input:
                cram_manifest = select_first([cram_manifest])
        }
    }

    Boolean do_cram_qc = select_first([ConvertCramManifestToCramArrays.passes_qc, defined(crams)]) # only do cram QC if manifest conversion passed QC
    
    # validations for array crams input
    if (do_cram_qc) {
        Array[File] cram_array = select_first([crams, ConvertCramManifestToCramArrays.crams])
        Array[File] cram_index_array = select_first([cram_indices, ConvertCramManifestToCramArrays.cram_indices])
        Array[String] sample_id_array = select_first([sample_ids, ConvertCramManifestToCramArrays.sample_ids])
        
        if (!defined(cram_index_array) || !defined(sample_id_array)) {
            Boolean no_cram_index_or_sample_id_passes_qc = false
            String no_cram_index_or_sample_id_message = "CRAM indices and sample IDs are required when CRAM files are provided. Please provide CRAM index files and a corresponding list of sample IDs."
        }

        call ValidateCramsAndIndices {
            input:
                crams = cram_array,
                cram_indices = cram_index_array,
                sample_ids = sample_id_array
        }
    }

    output {
        Boolean passes_qc = select_first([no_data_passes_qc, multiple_data_types_passes_qc, no_cram_index_or_sample_id_passes_qc, ValidateCramsAndIndices.passes_qc, ConvertCramManifestToCramArrays.passes_qc])
        String qc_messages = select_first([no_data_message, multiple_data_types_message, no_cram_index_or_sample_id_message, ValidateCramsAndIndices.qc_messages, ConvertCramManifestToCramArrays.qc_messages])
    }
}


task ConvertCramManifestToCramArrays {
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

        # convert the cram manifest into arrays of crams, cram indices, and sample ids
        head -n 1 ~{cram_manifest} > header.txt

        sample_id_col=$(cat header.txt | tr '\t' '\n' | grep -n "^sample_id$" | cut -d: -f1)
        cram_path_col=$(cat header.txt | tr '\t' '\n' | grep -n "^cram_path$" | cut -d: -f1)
        cram_index_col=$(cat header.txt | tr '\t' '\n' | grep -n "^cram_index_path$" | cut -d: -f1)

        if [ -z "${sample_id_col}" ] || [ -z "${cram_path_col}" ] || [ -z "${cram_index_col}" ]; then
            echo "Unable to determine column positions for sample_id, cram_path, or cram_index_path in the CRAM manifest." >> qc_messages.txt
            echo "false" > passes_qc.txt
            # create empty output files to ensure the task succeeds
            touch crams.txt cram_indices.txt sample_ids.txt
            exit 0
        fi

        tail -n +2 ~{cram_manifest} | cut -f${sample_id_col} > sample_ids.txt
        tail -n +2 ~{cram_manifest} | cut -f${cram_path_col} > crams.txt
        tail -n +2 ~{cram_manifest} | cut -f${cram_index_col} > cram_indices.txt

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
        Array[String] crams = read_lines("crams.txt")
        Array[String] cram_indices = read_lines("cram_indices.txt")
        Array[String] sample_ids = read_lines("sample_ids.txt")
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}

task ValidateCramsAndIndices {
    input {
        Array[String] crams
        Array[String] cram_indices
        Array[String] sample_ids

        Int max_cram_file_size_gb = 10

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = 10
    }

    Int num_crams = length(crams)
    Int num_cram_indices = length(cram_indices)
    Int num_sample_ids = length(sample_ids)

    command <<<
        # create empty qc messages file
        touch qc_messages.txt

        # write WDL arrays to files
        printf "%s\n" ~{sep=' ' sample_ids} > sample_ids_list.txt
        printf "%s\n" ~{sep=' ' crams} > crams_list.txt
        printf "%s\n" ~{sep=' ' cram_indices} > cram_indices_list.txt

        # validate that the number of CRAMs, CRAIs, and sample IDs match
        if [ ~{num_crams} -ne ~{num_cram_indices} ] || [ ~{num_crams} -ne ~{num_sample_ids} ]; then
            echo "Found different numbers of CRAMs (~{num_crams}), CRAIs (~{num_cram_indices}), and sample IDs (~{num_sample_ids})." >> qc_messages.txt
        else
            echo "Number of CRAMs, CRAIs, and sample IDs match: found ~{num_crams} of each."
        fi

        # validate that sample IDs are unique
        unique_sample_ids=$(cat sample_ids_list.txt | sort -u | wc -l)
        if [ $unique_sample_ids -ne ~{num_sample_ids} ]; then
            # find duplicate sample IDs
            duplicate_sample_ids=$(cat sample_ids_list.txt | sort | uniq -d | paste -sd, | sed 's/,/, /g')
            echo "Duplicate sample IDs found: ${duplicate_sample_ids}" >> qc_messages.txt
        else
            echo "Sample IDs are unique."
        fi

        # ensure all crams end with .cram and all cram indices end with .crai
        crams_with_wrong_extension=$(cat crams_list.txt | grep -vE "\.cram$" || true)
        if [ ! -z "${crams_with_wrong_extension}" ]; then
            echo "The following CRAM files do not have a .cram extension: ${crams_with_wrong_extension}" >> qc_messages.txt
        else
            echo "All CRAM files have the correct .cram extension."
        fi
        cram_indices_with_wrong_extension=$(cat cram_indices_list.txt | grep -vE "\.crai$" || true)
        if [ ! -z "${cram_indices_with_wrong_extension}" ]; then
            echo "The following CRAM index files do not have a .crai extension: ${cram_indices_with_wrong_extension}" >> qc_messages.txt
        else
            echo "All CRAM index files have the correct .crai extension."
        fi

        # validate that cram paths are unique
        unique_crams=$(cat crams_list.txt | sort -u | wc -l)
        if [ $unique_crams -ne ~{num_crams} ]; then
            # find duplicate CRAM paths
            duplicate_crams=$(cat crams_list.txt | sort | uniq -d | paste -sd, | sed 's/,/, /g')
            echo "Duplicate CRAM paths found: ${duplicate_crams}" >> qc_messages.txt
        else
            echo "CRAM paths are unique."
        fi

        # ensure that all CRAM files are less than the maximum file size allowed by the service (currently 10GB)
        # this also serves as an access check, which should already have been performed by the service
        crams_exceeding_max_size=$(cat crams_list.txt | while read cram; do
            file_size_bytes=$(gcloud storage ls -L "$cram" | grep "Content-Length:" | awk '{print $2}')
            file_size_gb=$((file_size_bytes / 1024 / 1024 / 1024))
            if [ $file_size_gb -gt ~{max_cram_file_size_gb} ]; then
                echo "$cram (${file_size_gb}GB)"
            fi
        done)
        if [ ! -z "${crams_exceeding_max_size}" ]; then
            echo "The following CRAM files exceed the maximum allowed file size of ~{max_cram_file_size_gb}GB: ${crams_exceeding_max_size}" >> qc_messages.txt
        else
            echo "All CRAM files are within the maximum allowed file size of ~{max_cram_file_size_gb}GB."
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
