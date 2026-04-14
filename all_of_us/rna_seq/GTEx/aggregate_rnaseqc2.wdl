version 1.0

task rnaseqc2_aggregate {
    input {
        File tpm_gcts_list         # newline-delimited list of TPM GCT file paths
        File count_gcts_list       # newline-delimited list of count GCT file paths
        File exon_count_gcts_list  # newline-delimited list of exon count GCT file paths
        File metrics_tsvs_list     # newline-delimited list of metrics TSV file paths
        String prefix
        File? insertsize_hists_list  # newline-delimited list of insert size histogram file paths
        String? flags

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command <<<
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Aggregating RNA-SeQC outputs")
        mkdir individual_outputs

        echo $(date +"[%b %d %H:%M:%S]") Staging TPM GCT files
        while read path; do
            gsutil -m cp "${path}" individual_outputs/
        done < ~{tpm_gcts_list}

        echo $(date +"[%b %d %H:%M:%S]") Staging count GCT files
        while read path; do
            gsutil -m cp "${path}" individual_outputs/
        done < ~{count_gcts_list}

        echo $(date +"[%b %d %H:%M:%S]") Staging exon count GCT files
        while read path; do
            gsutil -m cp "${path}" individual_outputs/
        done < ~{exon_count_gcts_list}

        echo $(date +"[%b %d %H:%M:%S]") Staging metrics TSV files
        while read path; do
            gsutil -m cp "${path}" individual_outputs/
        done < ~{metrics_tsvs_list}

        if [ -f ~{insertsize_hists_list} ]; then
            echo $(date +"[%b %d %H:%M:%S]") Staging insert size histogram files
            while read path; do
                gsutil -m cp "${path}" individual_outputs/
            done < ~{insertsize_hists_list}
        fi

        touch ~{prefix}.insert_size_hists.txt.gz
        python3 -m rnaseqc aggregate \
            -o . \
            individual_outputs \
            ~{prefix} \
            ~{flags}
        echo $(date +"[%b %d %H:%M:%S] done")
    >>>

    output {
        File metrics = "~{prefix}.metrics.txt.gz"
        File insert_size_hists = "~{prefix}.insert_size_hists.txt.gz"
        File tpm_gct = glob("~{prefix}.gene_tpm.*")[0]
        File count_gct = glob("~{prefix}.gene_reads.*")[0]
        File exon_count_gct = glob("~{prefix}.exon_reads.*")[0]
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


workflow rnaseqc2_aggregate_workflow {
    input {
        File tpm_gcts_list
        File count_gcts_list
        File exon_count_gcts_list
        File metrics_tsvs_list
        String prefix
        File? insertsize_hists_list
        String? flags
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }
	String pipeline_version = "aou_9.0.0"

    call rnaseqc2_aggregate {
        input:
            tpm_gcts_list = tpm_gcts_list,
            count_gcts_list = count_gcts_list,
            exon_count_gcts_list = exon_count_gcts_list,
            metrics_tsvs_list = metrics_tsvs_list,
            prefix = prefix,
            insertsize_hists_list = insertsize_hists_list,
            flags = flags,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            num_preempt = num_preempt
	}
 

    output {
        File metrics = rnaseqc2_aggregate.metrics
        File insert_size_hists = rnaseqc2_aggregate.insert_size_hists
        File tpm_gct = rnaseqc2_aggregate.tpm_gct
        File count_gct = rnaseqc2_aggregate.count_gct
        File exon_count_gct = rnaseqc2_aggregate.exon_count_gct
    }
}