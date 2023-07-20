version 1.0

workflow NewAndImprovedsn3MC {

  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
  }

  call sort_r1_and_r2 {
    input:
      r1 = fastq_input_read1,
      r2 = fastq_input_read2
  }

  call trim {
    input:
     r1_sorted = sort_r1_and_r2.r1_sorted_fq,
     r2_sorted = sort_r1_and_r2.r2_sorted_fq
  }

  call  hisat_3n_pair_end_mapping_dna_mode {
    input:
     r1_trimmed = trim.r1_trimmed_fq,
     r2_trimmed = trim.r2_trimmed_fq
  }

  call separate_unmapped_reads {
    input:
      hisat3n_bam = hisat_3n_pair_end_mapping_dna_mode.hisat3n_bam
  }

  call split_unmapped_reads {
    input:
      unmapped_fastq = separate_unmapped_reads.unmapped_fastq
  }

  call hisat_3n_single_end_r1_mapping_dna_mode {
    input:
      split_r1 = split_unmapped_reads.split_r1_fq
  }

  call hisat_3n_single_end_r2_mapping_dna_mode {
    input:
      split_r2 = split_unmapped_reads.split_r2_fq
  }

  call merge_and_sort_split_reads_by_name {
    input:
      r1_hisat3n_bam = hisat_3n_single_end_r1_mapping_dna_mode.r1_hisat3n_bam,
      r2_hisat3n_bam = hisat_3n_single_end_r2_mapping_dna_mode.r2_hisat3n_bam
  }

  call remove_overlap_read_parts {
    input:
      bam = merge_and_sort_split_reads_by_name.bam
  }

  call merge_original_and_split_bam {
    input:
      bam = separate_unmapped_reads.unique_bam,
      split_bam = remove_overlap_read_parts.remove_overlap_bam
  }

  call sort_all_reads_by_name {
    input:
      bam = merge_original_and_split_bam.bam
  }

  call call_chromatin_contacts {
    input:
      bam = sort_all_reads_by_name.bam
  }

  call sort_bam {
    input:
      bam = sort_all_reads_by_name.bam
  }

  call dedup_unique_bam {
    input:
      bam = sort_bam.bam
  }

  call index_unique_bam_dna_reads {
    input:
      bam = dedup_unique_bam.dedup_bam
  }

  call unique_reads_allc {
    input:
      bam = dedup_unique_bam.dedup_bam,
      bai = index_unique_bam_dna_reads.bai
  }

  call unique_reads_cgn_extraction {
    input:
      allc = unique_reads_allc.allc,
      tbi = unique_reads_allc.tbi
  }

  call summary {
    input:
      trimmed_stats = trim.trim_stats,
      hisat3n_stats = hisat_3n_pair_end_mapping_dna_mode.hisat3n_stats,
      r1_hisat3n_stats = hisat_3n_single_end_r1_mapping_dna_mode.r1_hisat3n_stats,
      r2_hisat3n_stats = hisat_3n_single_end_r2_mapping_dna_mode.r2_hisat3n_stats,
      dedup_stats = dedup_unique_bam.dedup_stats,
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


task sort_r1_and_r2 {
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
        File r1_sorted_fq = ""
        File r2_sorted_fq = ""
    }
}

task trim {
    input {
    File r1_sorted
    File r2_sorted
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

task hisat_3n_single_end_r1_mapping_dna_mode {
    input {
           File split_r1
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
    }
}

task hisat_3n_single_end_r2_mapping_dna_mode {
    input {
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
        File r2_hisat3n_bam = ""
        File r2_hisat3n_stats = ""
    }
}

task merge_and_sort_split_reads_by_name {
    input {
        File r1_hisat3n_bam
        File r2_hisat3n_bam
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
        File bam = ""
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

task merge_original_and_split_bam {
    input {
        File bam
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
        File bam = ""
    }
}

task sort_all_reads_by_name {
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
        File bam = ""
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

task sort_bam {
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
        File bam = ""
    }
}

task dedup_unique_bam {
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
    }
}

task index_unique_bam_dna_reads {
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
        File bai = ""
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