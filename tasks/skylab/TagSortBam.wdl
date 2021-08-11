version 1.0

task CellSortBam {
  input {
    File bam_input
    File original_gtf

    # runtime values
    String docker = "quay.io/kishorikonwar/secondary-analysis-python3-scientific:sctools-optimized8"
    Int machine_mem_mb = 8000
    Int cpu = 8
    Int disk = ceil(size(bam_input, "Gi") * 4)
    Int preemptible = 0
  }

  meta {
    description: "Sort bam_input by cell, then molecule, then gene."
  }

  parameter_meta {
    bam_input: "Input bam file containing reads marked with tags for cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  
  command {
    set -e
    mkdir temp
    gunzip -c ~{original_gtf} > annotation.gtf

    TagSort --bam-input ~{bam_input} \
    --gtf-file annotation.gtf \
    --metric-output cell-metrics.csv  --compute-metric \
    --metric-type cell \
    --barcode-tag CB \
    --umi-tag UB \
    --gene-tag GX \
    --temp-folder temp \
    --alignments-per-thread 1000000 \
    --nthreads 6

    gzip cell-metrics.csv
  }


  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File cell_metrics = "cell-metrics.csv.gz"
  }
}

task GeneSortBam {
  input {
    File bam_input

    # runtime values
    String docker = "quay.io/kishorikonwar/secondary-analysis-python3-scientific:sctools-optimized8"
    Int machine_mem_mb = 8000
    Int cpu = 8
    Int disk = ceil(size(bam_input, "Gi") * 4) 
    Int preemptible = 0
  }
  

  meta {
    description: "Sort bam_input by gene, then cell, then molecule."
  }

  parameter_meta {
    bam_input: "Input bam file containing reads marked with tags for cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e
    mkdir temp

    TagSort --bam-input ~{bam_input} \
    --metric-output gene-metrics.csv  --compute-metric \
    --metric-type gene \
    --gene-tag GX \
    --barcode-tag CB \
    --umi-tag UB \
    --temp-folder temp \
    --alignments-per-thread 1000000 \
    --nthreads 6

    gzip gene-metrics.csv

  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD" #TODO: make it an input
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File gene_metrics = "gene-metrics.csv.gz"
  }
}
