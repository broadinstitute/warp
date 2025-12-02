version 1.0
task leafcutter_cluster {
	input {

		File junc_files_list
		File exon_list
		File genes_gtf
		String prefix
		File sample_participant_lookup

		Int? min_clu_reads
		Float? min_clu_ratio
		Int? max_intron_len
		Int? num_pcs

		Int memory
		Int disk_space
		Int num_threads
		Int num_preempt
	}


    command <<<
        export PATH=$PATH:/root/google-cloud-sdk/bin
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S]") Starting leafcutter
        # Tune gsutil parallelism (adjust to your VM)
        export GSUTIL_PARALLEL_PROCESS_COUNT=32
        export GSUTIL_PARALLEL_THREAD_COUNT=8

        mkdir -p junc_inputs

        echo $(date +"[%b %d %H:%M:%S]") "[leafcutter] parallel copy begin"
        # Single invocation, reads entire list from STDIN
        gsutil -m cp -I junc_inputs/ < ~{junc_files_list}
        echo $(date +"[%b %d %H:%M:%S]") "[leafcutter] parallel copy done"

        # Build the list for the python script
        junc_files=(junc_inputs/*)
        printf "%s\n" "${junc_files[@]}" > junc_file_list.txt

        # mkdir junc_inputs
        # while read path; do
        #     gsutil -m cp "${path}" junc_inputs/
        # done < ~{junc_files_list}

        # junc_files=(junc_inputs/*)
        printf "%s\n" ${junc_files[@]} > junc_file_list.txt		
        python3 /src/cluster_prepare_fastqtl.py \
            junc_file_list.txt \
            ~{exon_list} \
            ~{genes_gtf} \
            ~{prefix} \
            ~{sample_participant_lookup} \
            ~{"--min_clu_reads " + min_clu_reads} \
            ~{"--min_clu_ratio " + min_clu_ratio} \
            ~{"--max_intron_len " + max_intron_len} \
            ~{"--num_pcs " + num_pcs} 
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/leafcutter:1.0.0"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File counts="~{prefix}_perind.counts.gz"
        File counts_numers="~{prefix}_perind_numers.counts.gz"
        File clusters_pooled="~{prefix}_pooled.gz"
        File clusters_refined="~{prefix}_refined.gz"
        File phenotype_groups="~{prefix}.leafcutter.phenotype_groups.txt"
        File leafcutter_bed_parquet="~{prefix}.leafcutter.bed.parquet"
        File leafcutter_bed="~{prefix}.leafcutter.bed.gz"
        File leafcutter_bed_index="~{prefix}.leafcutter.bed.gz.tbi"
        File leafcutter_pcs="~{prefix}.leafcutter.PCs.txt"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow leafcutter_cluster_workflow {
	input {
		File junc_files_list
		File exon_list
		File genes_gtf
		String prefix
		File sample_participant_lookup

		Int? min_clu_reads
		Float? min_clu_ratio
		Int? max_intron_len
		Int? num_pcs

		Int memory
		Int disk_space
		Int num_threads
		Int num_preempt
	}
	String pipeline_version = "aou_9.0.0"
    call leafcutter_cluster {
		input:
			junc_files_list = junc_files_list,
			exon_list = exon_list,
			genes_gtf = genes_gtf,
			prefix = prefix,
			sample_participant_lookup = sample_participant_lookup,
			min_clu_reads = min_clu_reads,
			min_clu_ratio = min_clu_ratio,
			max_intron_len = max_intron_len,
			num_pcs = num_pcs,
			memory = memory,
			disk_space = disk_space,
			num_threads = num_threads,
			num_preempt = num_preempt
	}
	output {
		File counts_out = leafcutter_cluster.counts
		File counts_numers_out = leafcutter_cluster.counts_numers
		File clusters_pooled_out = leafcutter_cluster.clusters_pooled
		File clusters_refined_out = leafcutter_cluster.clusters_refined
		File phenotype_groups_out = leafcutter_cluster.phenotype_groups
		File bed_parquet_out = leafcutter_cluster.leafcutter_bed_parquet
		File bed_out = leafcutter_cluster.leafcutter_bed
		File bed_index_out = leafcutter_cluster.leafcutter_bed_index
		File pcs_out = leafcutter_cluster.leafcutter_pcs
		String leafcutter_pipeline_version = pipeline_version
	}
}