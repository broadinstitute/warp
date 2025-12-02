version 1.0
task rsem_aggregate_results {
    input {

        File rsem_isoforms_list  # newline-delimited list of isoform result paths
        File rsem_genes_list     # newline-delimited list of gene result paths
        String prefix

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
        # Captures version used for AoU processing
        String pipeline_version = "aou_9.0.0"

    }
    command <<<
        echo $(date +"[%b %d %H:%M:%S]") Starting transcript-level aggregation

        mkdir isoform_inputs
        while read path; do
            gsutil -m cp "${path}" isoform_inputs/
        done < ~{rsem_isoforms_list}

        isoform_files=(isoform_inputs/*)
        printf "%s\n" ${isoform_files[@]} > isoform_file_list.txt
        python3 /src/aggregate_rsem_results.py isoform_file_list.txt TPM IsoPct expected_count ~{prefix}
        
        echo $(date +"[%b %d %H:%M:%S]") Starting gene-level aggregation

        mkdir gene_inputs
        while read path; do
            gsutil -m cp "${path}" gene_inputs/
        done < ~{rsem_genes_list}

        gene_files=(gene_inputs/*)
        printf "%s\n" ${gene_files[@]} > gene_file_list.txt
        python3 /src/aggregate_rsem_results.py gene_file_list.txt TPM expected_count ~{prefix}
    >>>

    output {
        File transcripts_tpm="${prefix}.rsem_transcripts_tpm.txt.gz"
        File transcripts_isopct="${prefix}.rsem_transcripts_isopct.txt.gz"
        File transcripts_expected_count="${prefix}.rsem_transcripts_expected_count.txt.gz"
        File genes_tpm="${prefix}.rsem_genes_tpm.txt.gz"
        File genes_expected_count="${prefix}.rsem_genes_expected_count.txt.gz"
        String pipeline_version_out = "~{pipeline_version}"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow aggregate_rsem_results {
    call rsem_aggregate_results
}