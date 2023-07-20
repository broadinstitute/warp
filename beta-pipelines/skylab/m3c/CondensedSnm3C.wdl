version 1.0

workflow NewAndImprovedsn3MC {

    input {
        Array[File] fastq_input_read1
        Array[File] fastq_input_read2
    }

    call sort_and_trim_r1_and_r2 {
        input:
            r1 = fastq_input_read1,
            r2 = fastq_input_read2
    }

    call hisat_3n_pair_end_mapping_dna_mode {
        input:
            r1_trimmed = sort_and_trim_r1_and_r2.r1_trimmed_fq,
            r2_trimmed = sort_and_trim_r1_and_r2.r2_trimmed_fq
    }

    call separate_unmapped_reads {
        input:
            hisat3n_bam = hisat_3n_pair_end_mapping_dna_mode.hisat3n_bam
    }

    call split_unmapped_reads {
        input:
            unmapped_fastq = separate_unmapped_reads.unmapped_fastq
    }

    call hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name {
        input:
            split_r1 = split_unmapped_reads.split_r1_fq,
            split_r2 = split_unmapped_reads.split_r2_fq
    }

    call remove_overlap_read_parts {
        input:
            bam = hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.merge_sorted_bam
    }

    call merge_original_and_split_bam_and_sort_all_reads_by_name_and_position {
        input:
            unique_bam = separate_unmapped_reads.unique_bam,
            split_bam = remove_overlap_read_parts.remove_overlap_bam
    }

    call call_chromatin_contacts {
        input:
            bam = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position.name_sorted_bam
    }

    call dedup_unique_bam_and_index_unique_bam {
        input:
            bam = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position.position_sorted_bam
    }

    call unique_reads_allc {
        input:
            bam = dedup_unique_bam_and_index_unique_bam.dedup_bam,
            bai = dedup_unique_bam_and_index_unique_bam.dedup_bai
    }

    call unique_reads_cgn_extraction {
        input:
            allc = unique_reads_allc.allc,
            tbi = unique_reads_allc.tbi
    }

    call summary {
        input:
            trimmed_stats = sort_and_trim_r1_and_r2.trim_stats,
            hisat3n_stats = hisat_3n_pair_end_mapping_dna_mode.hisat3n_stats,
            r1_hisat3n_stats = hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.r1_hisat3n_stats,
            r2_hisat3n_stats = hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.r2_hisat3n_stats,
            dedup_stats = dedup_unique_bam_and_index_unique_bam.dedup_stats,
            chromatin_contact_stats = call_chromatin_contacts.chromatin_contact_stats,
            allc_uniq_reads_stats = unique_reads_allc.allc_uniq_reads_stats,
            unique_reads_cgn_extraction_tbi = unique_reads_cgn_extraction.unique_reads_cgn_extraction_tbi
    }

    output {
        File MappingSummary = summary.MappingSummary
        #File allcFiles
        #File allc_CGNFiles
        #File bamFiles
        #File detail_statsFiles =
        #File hicFiles
    }
}

task sort_and_trim_r1_and_r2 {
    input {
        File r1
        File r2
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File r1_trimmed_fq = ""
        File r2_trimmed_fq = ""
        File trim_stats = ""
    }
}

task hisat_3n_pair_end_mapping_dna_mode{
    input {
        File r1_trimmed
        File r2_trimmed
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File hisat3n_bam = ""
        File hisat3n_stats = ""
    }
}

task separate_unmapped_reads {
    input {
        File hisat3n_bam
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File unique_bam = ""
        File multi_bam = ""
        File unmapped_fastq = ""
    }
}

task split_unmapped_reads {
    input {
        File unmapped_fastq
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File split_r1_fq = ""
        File split_r2_fq = ""
    }
}

task hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name {
    input {
        File split_r1
        File split_r2
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File r1_hisat3n_bam = ""
        File r1_hisat3n_stats = ""
        File r2_hisat3n_bam = ""
        File r2_hisat3n_stats = ""
        File merge_sorted_bam = ""
    }
}

task remove_overlap_read_parts {
    input {
        File bam
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File remove_overlap_bam = ""
    }
}

task merge_original_and_split_bam_and_sort_all_reads_by_name_and_position {
    input {
        File unique_bam
        File split_bam
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File name_sorted_bam = ""
        File position_sorted_bam = ""
    }
}

task call_chromatin_contacts {
    input {
        File bam
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File chromatin_contact_stats = ""
    }
}

task dedup_unique_bam_and_index_unique_bam {
    input {
        File bam
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File dedup_bam = ""
        File dedup_stats = ""
        File dedup_bai = ""
    }
}

task unique_reads_allc {
    input {
        File bam
        File bai
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File allc = ""
        File tbi = ""
        File allc_uniq_reads_stats = ""
    }
}

task unique_reads_cgn_extraction {
    input {
        File allc
        File tbi
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File unique_reads_cgn_extraction_allc = ""
        File unique_reads_cgn_extraction_tbi = ""
    }
}

task summary {
    input {
        File trimmed_stats
        File hisat3n_stats
        File r1_hisat3n_stats
        File r2_hisat3n_stats
        File dedup_stats
        File chromatin_contact_stats
        File allc_uniq_reads_stats
        File unique_reads_cgn_extraction_tbi
    }
    command <<<
    >>>
    runtime {
        docker: "fill_in"
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File mapping_summary = "MappingSummary.csv.gz"
    }
}