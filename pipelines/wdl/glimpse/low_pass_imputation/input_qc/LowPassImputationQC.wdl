version 1.0

workflow InputQC {
    input {
        String pipeline_version = "0.0.1"

        # List of files, one per line
        File reference_chunks
        File sites_vcf
        File sites_table
        File sites_table_index

        Array[String] contigs

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

    if (defined(crams)) {
        if (!defined(cram_indices) || !defined(sample_ids)) {
            Boolean no_cram_index_or_sample_id_passes_qc = false
            String no_cram_index_or_sample_id_message = "CRAM indices and sample IDs are required when CRAM files are provided. Please provide cram index files and a list of sample IDs corresponding to the CRAM files."
        }
        # call ValidateCramsAndIndices {
        #     input:
        #         crams = crams,
        #         cram_indices = cram_indices,
        #         sample_ids = sample_ids,
        #         cram_manifest = cram_manifest
        # }
    }

    if (defined(cram_manifest)) {
        call ValidateCramManifest {
            input:
                cram_manifest = select_first([cram_manifest])
        }
    }

    output {
        Boolean passes_qc = select_first([no_data_passes_qc, multiple_data_types_passes_qc, no_cram_index_or_sample_id_passes_qc])
        String qc_messages = select_first([no_data_message, multiple_data_types_message, no_cram_index_or_sample_id_message])
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

        # validate that the manifest has the required columns
        required_columns="sample_id\tcram_path\tcram_index_path"
        header=$(head -n 1 ${cram_manifest})
        if [[ "${header}" != "${required_columns}" ]]; then
            echo "CRAM manifest is missing required columns. Expected header: ${required_columns}. Actual header: ${header}" >> qc_messages.txt
        else
            echo "CRAM manifest has the required columns."
        fi

        # TODO add more checks

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
