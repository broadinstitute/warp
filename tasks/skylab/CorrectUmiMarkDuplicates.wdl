task SortAndCorrectUmiMarkDuplicates {
  File bam_input

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
  Int machine_mem_mb = 7500
  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500
  Int cpu = 1
  # mark duplicates swaps a large amount of data to disk, has high disk requirements.
  Int disk = ceil(size(bam_input, "G") * 6) + 50
  Int preemptible = 3

  meta {
    description: "Iterates over reads in a bam file, marking duplicate reads by setting bit 0x400 in the SAM flag and correcting any barcodes with errors by storing the original barcode in tag UR."
  }

  parameter_meta {
    bam_input: "Aligned bam file containing reads that have been marked by cell barcodes (CB) and molecular barcodes (UB)"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar SortSam \
      SORT_ORDER=coordinate \
      I=${bam_input} \
      O=sorted.bam

    # recover disk space
    rm ${bam_input}

    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar UmiAwareMarkDuplicatesWithMateCigar \
      MAX_EDIT_DISTANCE_TO_JOIN=1 \
      UMI_METRICS_FILE=umi_metrics.txt \
      UMI_TAG_NAME=UR \
      ASSIGNED_UMI_TAG=UB \
      BARCODE_TAG=CB \
      METRICS_FILE=duplicate_metrics.txt \
      OUTPUT=duplicates_marked.bam \
      INPUT=sorted.bam
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "duplicates_marked.bam"
    File umi_metrics = "umi_metrics.txt"
    File duplicate_metrics = "duplicate_metrics.txt"
  }
} 
