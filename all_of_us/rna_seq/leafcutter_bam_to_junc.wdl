version 1.0
task leafcutter_bam_to_junc {
	input {
		File bam_file
		String sample_id
		Int? strand_specificity = 0

		Int memory
		Int disk_space
		Int num_threads
		Int num_preempt
	}


    command <<<
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Extracting junctions for sample ~{sample_id}")
        # select uniquely mapped reads that pass WASP filters
        filtered_bam=~{bam_file}.filtered.bam
        samtools view -h -q 255 ~{bam_file} | grep -v "vW:i:[2-7]" | samtools view -b > $filtered_bam
        samtools index $filtered_bam
        regtools junctions extract -a 8 -m 50 -M 500000 -s ~{strand_specificity} $filtered_bam | gzip -c > ~{sample_id}.regtools_junc.txt.gz
        echo $(date +"[%b %d %H:%M:%S] Done")
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/leafcutter:latest"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File junc_file="~{sample_id}.regtools_junc.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow leafcutter_bam_to_junc_workflow {
    input {
		File bam_file
		String sample_id
		Int? strand_specificity
        Int? user_defined_threads = 1  # optional override
        Int? user_defined_preempt = 1  # optional override
	}
	String pipeline_version = "aou_9.0.0"
	   # Get size of BAM file in GiB (1 GiB = 1024^3 bytes)
    Float bam_file_size_gb = size(bam_file, "GB")

    # Compute runtime settings based on file size
    Int memory = ceil(bam_file_size_gb * 2.0) + 1
    Int disk_space = ceil(bam_file_size_gb * 3.0) + 5
    Int num_threads = select_first([user_defined_threads, 1])
    Int num_preempt = select_first([user_defined_preempt, 1])
	
	call leafcutter_bam_to_junc {
		input:
			bam_file = bam_file,
			sample_id = sample_id,
			strand_specificity = strand_specificity,

			memory = memory,
			disk_space = disk_space,
			num_threads = num_threads,
			num_preempt = num_preempt
	}
	output {
		File junc_file_final = leafcutter_bam_to_junc.junc_file
		String pipeline_version_out = pipeline_version
	}
}
