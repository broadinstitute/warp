version 1.0

task SortBam {
    input {
        File bam_input
        String sort_order = "coordinate"

        # runtime values
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        Int machine_mem_mb = 8250
        Int machine_overhead_mb = 500
        Int cpu = 1
        Int preemptible = 3
    }

    Int command_mem_mb = machine_mem_mb - machine_overhead_mb
    Int disk = ceil(size(bam_input, "Gi") * 6) + 200

    meta {
        description: "Sorts bam"
    }

    command {
        set -e

        java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar SortSam \
              I=${bam_input} \
              O=sorted.bam \
              SORT_ORDER=${sort_order}
    }

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File bam_output = "sorted.bam"
    }
}

task SortBamAndIndex {
    input {
        File bam_input
        String sort_order = "coordinate"

        # runtime values
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        Int machine_mem_mb = 8250
        Int machine_overhead_mb = 500
        Int command_mem_mb = machine_mem_mb - machine_overhead_mb
        Int cpu = 1
        Int disk = ceil(size(bam_input, "Gi") * 6) + 200
        Int preemptible = 3
    }

    meta {
        description: "Sorts bam by genomic position"
    }

    command {
        set -e
        java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar SortSam \
              I=${bam_input} \
              O=sorted.bam \
              SORT_ORDER=${sort_order}
         java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar BuildBamIndex \
              I=sorted.bam \
              O=sorted.bai
    }

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File bam_output = "sorted.bam"
        File bam_index = "sorted.bai"
    }
}

task CollectMultipleMetrics {
  input {
    File aligned_bam
    File genome_ref_fasta
    String output_basename

    # runtime values
    String docker ="us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    Int machine_mem_mb = 8250
    # give the command 1 GiB of overhead
    Int command_mem_mb = machine_mem_mb - 1000
    Int cpu = 1
    # use provided disk number or dynamically size on our own, with 200GiB of additional disk
    Int disk = ceil(size(aligned_bam, "GiB") + size(genome_ref_fasta, "GiB") + 200)
    Int preemptible = 3
  }

  meta {
    description: "This Picard task will collect multiple QC metrics, such as CollectAlignmentSummaryMetrics and CollectInsertSizeMetrics."
  }

  parameter_meta {
    aligned_bam: "input aligned bam"
    genome_ref_fasta: "genome reference fasta"
    output_basename: "basename used for output files"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_basename}" \
      FILE_EXTENSION=".txt" \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectGcBiasMetrics \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=QualityScoreDistribution \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=CollectQualityYieldMetrics \
      REFERENCE_SEQUENCE="${genome_ref_fasta}" \
      ASSUME_SORTED=true
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File alignment_summary_metrics = "${output_basename}.alignment_summary_metrics.txt"
    File base_call_dist_metrics = "${output_basename}.base_distribution_by_cycle_metrics.txt"
    File base_call_pdf = "${output_basename}.base_distribution_by_cycle.pdf"
    File gc_bias_detail_metrics = "${output_basename}.gc_bias.detail_metrics.txt"
    File gc_bias_dist_pdf = "${output_basename}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_basename}.gc_bias.summary_metrics.txt"
    Array[File] insert_size_hist = glob("${output_basename}.insert_size_histogram.pdf")
    Array[File?] insert_size_metrics = glob("${output_basename}.insert_size_metrics.txt")
    File quality_distribution_metrics = "${output_basename}.quality_distribution_metrics.txt"
    File quality_distribution_dist_pdf = "${output_basename}.quality_distribution.pdf"
    File quality_by_cycle_metrics = "${output_basename}.quality_by_cycle_metrics.txt"
    File quality_by_cycle_pdf = "${output_basename}.quality_by_cycle.pdf"
    File pre_adapter_details_metrics = "${output_basename}.pre_adapter_detail_metrics.txt"
    File pre_adapter_summary_metrics = "${output_basename}.pre_adapter_summary_metrics.txt"
    File bait_bias_detail_metrics = "${output_basename}.bait_bias_detail_metrics.txt"
    File bait_bias_summary_metrics = "${output_basename}.bait_bias_summary_metrics.txt"
    File error_summary_metrics = "${output_basename}.error_summary_metrics.txt"
  }
}
task CollectMultipleMetricsMultiSample {
    input {
        Array[File] aligned_bam_inputs
        File genome_ref_fasta
        Array[String] input_ids

        # runtime values
        String docker ="us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        Int machine_mem_mb = 8250
        # give the command 1 GiB of overhead
        Int command_mem_mb = machine_mem_mb - 1000
        Int cpu = 4
        # use provided disk number or dynamically size on our own, with 200GiB of additional disk
        Int disk = ceil(size(aligned_bam_inputs, "GiB") + size(genome_ref_fasta, "GiB") + 50)
        Int preemptible = 3
    }

    meta {
        description: "This Picard task will collect multiple QC metrics, such as CollectAlignmentSummaryMetrics and CollectInsertSizeMetrics."
    }

    parameter_meta {
        aligned_bam_inputs: "Array of input aligned bam files"
        genome_ref_fasta: "genome reference fasta"
        input_ids: "basename used for output files"
        docker: "(optional) the docker image containing the runtime environment for this task"
        machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command <<<
        set -e

        declare -a bam_files=(~{sep=' ' aligned_bam_inputs})
        declare -a output_prefix=(~{sep=' ' input_ids})
        for (( i=0; i<${#bam_files[@]}; ++i));
        do
            output_basename=${output_prefix[$i]}
            java -Xmx"~{command_mem_mb}"m \
            -jar /usr/picard/picard.jar CollectMultipleMetrics \
            --VALIDATION_STRINGENCY SILENT \
            --METRIC_ACCUMULATION_LEVEL ALL_READS \
            --INPUT "${bam_files[$i]}" \
            --OUTPUT "${output_basename}" \
            --FILE_EXTENSION ".txt" \
            --PROGRAM null \
            --PROGRAM CollectAlignmentSummaryMetrics  \
            --PROGRAM CollectGcBiasMetrics \
            --REFERENCE_SEQUENCE "~{genome_ref_fasta}" \
            --ASSUME_SORTED true
        done;
    >>>

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        cpu: cpu
        preemptible: preemptible
    }

    output {
        Array[File] alignment_summary_metrics = glob("*.alignment_summary_metrics.txt")
        Array[File] gc_bias_summary_metrics = glob("*.gc_bias.summary_metrics.txt")

    }
}
task CollectRnaMetrics {
  input {
    File aligned_bam
    File ref_flat
    File rrna_intervals
    String output_basename
    String stranded

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    Int machine_mem_mb = 3850
    # give the command 500 MiB of overhead
    Int command_mem_mb = machine_mem_mb - 500
    Int cpu = 1
    # use provided disk number or dynamically size on our own, with 200GiB of additional disk
    Int disk = ceil(size(aligned_bam, "GiB") + size(ref_flat, "GiB") + size(rrna_intervals, "GiB") + 200)
    Int preemptible = 3
  }
  

  meta {
    description: "This Picard task will collect RnaSeqMetrics."
  }

  parameter_meta {
    aligned_bam: "input aligned file"
    ref_flat: "reference flat file"
    rrna_intervals: "ribosomal intervals"
    output_basename: "basename used for output files"
    stranded: "whether or not to use strand specificity"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  
  command {
    set -e
    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_basename}.rna_metrics.txt" \
      REF_FLAT="${ref_flat}" \
      RIBOSOMAL_INTERVALS="${rrna_intervals}" \
      STRAND_SPECIFICITY=${stranded} \
      CHART_OUTPUT="${output_basename}.rna.coverage.pdf"
    touch "${output_basename}.rna.coverage.pdf"
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File rna_metrics = "${output_basename}.rna_metrics.txt"
    File rna_coverage_pdf = "${output_basename}.rna.coverage.pdf"
  }
}

# Here we use "-XX:ParallelGCThreads=2" to run MarkDuplication on multiple threads 
task CollectDuplicationMetrics {
  input {
    File aligned_bam
    String output_basename

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    Int machine_mem_mb = 32768
    # give the command 1 GiB of overhead
    Int command_mem_mb = machine_mem_mb - 1000
    Int cpu = 2
    # use provided disk number or dynamically size on our own, with 200GiB of additional disk
    Int disk = ceil(size(aligned_bam, "GiB") + 200)
    Int preemptible = 3
  }
  

  meta {
    description: "This Picard task will collect alignment DuplicationMetrics."
  }

  parameter_meta {
    aligned_bam: "input aligned bam"
    output_basename: "basename used for output files"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  
  command {
    java -Xmx${command_mem_mb}m -XX:ParallelGCThreads=${cpu}  -jar /usr/picard/picard.jar  MarkDuplicates \
       VALIDATION_STRINGENCY=SILENT  \
       INPUT=${aligned_bam} \
       OUTPUT="${output_basename}.MarkDuplicated.bam" \
       ASSUME_SORTED=true \
       METRICS_FILE="${output_basename}.duplicate_metrics.txt" \
       REMOVE_DUPLICATES=false
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File dedup_metrics = "${output_basename}.duplicate_metrics.txt"
  }
}

# Here we use "-XX:ParallelGCThreads=2" to run MarkDuplication on multiple threads 
task RemoveDuplicatesFromBam {
  input {
    Array[File] aligned_bam_inputs
    Array[String] input_ids

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    Int machine_mem_mb = 32768
    # give the command 1 GiB of overhead
    Int command_mem_mb = machine_mem_mb - 1000
    Int cpu = 2
    # use provided disk number or dynamically size on our own, with 200GiB of additional disk
    Int disk = ceil(size(aligned_bam_inputs, "GiB") * 2.5) + 10
    Int preemptible = 3
  }
  

  meta {
    description: "This Picard task will collect alignment DuplicationMetrics."
  }

  parameter_meta {
    aligned_bam_inputs: "input aligned bam"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -e

    declare -a bam_files=(~{sep=' ' aligned_bam_inputs})
    declare -a output_prefix=(~{sep=' ' input_ids})
    for (( i=0; i<${#bam_files[@]}; ++i));
    do
      java -Xmx"~{command_mem_mb}"m -XX:ParallelGCThreads=~{cpu} -jar /usr/picard/picard.jar  MarkDuplicates \
       -VALIDATION_STRINGENCY SILENT  \
       -INPUT "${bam_files[$i]}" \
       -OUTPUT "${output_prefix[$i]}.aligned_bam.DuplicatesRemoved.bam" \
       -ASSUME_SORTED true \
       -METRICS_FILE "${output_prefix[$i]}.aligned_bam.duplicate_metrics.txt" \
       -REMOVE_DUPLICATES true;

    java -Xmx"~{command_mem_mb}"m -XX:ParallelGCThreads=~{cpu} -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
       -I "${output_prefix[$i]}.aligned_bam.DuplicatesRemoved.bam" \
       -O "${output_prefix[$i]}.aligned_bam.DuplicatesRemoved.ReadgroupAdded.bam" \
       -RGID 4 \
       -RGLB lib1 \
       -RGPL ILLUMINA \
       -RGPU unit1 \
       -RGSM 20
      done;
  >>>
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    Array[File] dedup_metrics = glob("*.aligned_bam.duplicate_metrics.txt")
    Array[File] output_bam = glob("*.aligned_bam.DuplicatesRemoved.ReadgroupAdded.bam")
  }
}

