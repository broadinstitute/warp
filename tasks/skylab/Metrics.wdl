version 1.0

task CalculateCellMetrics {
  input {
    File bam_input
    File original_gtf

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.4.0"
    Int machine_mem_mb = 8000
    Int cpu = 4
    Int disk = ceil(size(bam_input, "Gi") * 4) + ceil((size(original_gtf, "Gi") * 3)) 
    Int preemptible = 3
  }

  meta {
    description: "Sort bam_input by cell, then molecule, then gene."
  }

  parameter_meta {
    bam_input: "Input bam file containing reads marked with tags for cell barcodes (CB), molecule barcodes (UB) and gene ids (GX)"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  
  command {
    set -e
    
    # create the tmp folder for disk sorting
    mkdir temp

    # if GTF file in compressed then uncompress
    if [[ ~{original_gtf} =~ \.gz$ ]]
    then
        gunzip -c ~{original_gtf} > annotation.gtf
    else
        mv  ~{original_gtf}  annotation.gtf
    fi

    TagSort --bam-input ~{bam_input} \
    --gtf-file annotation.gtf \
    --metric-output cell-metrics.csv \
    --compute-metric \
    --metric-type cell \
    --barcode-tag CB \
    --umi-tag UB \
    --gene-tag GX \
    --temp-folder temp \
    --alignments-per-thread 1000000 \
    --nthreads ${cpu}

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

task CalculateGeneMetrics {
  input {
    File bam_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.4.0"
    Int machine_mem_mb = 8000
    Int cpu = 4
    Int disk = ceil(size(bam_input, "Gi") * 4) 
    Int preemptible = 3
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
    --metric-output gene-metrics.csv \
    --compute-metric \
    --metric-type gene \
    --gene-tag GX \
    --barcode-tag CB \
    --umi-tag UB \
    --temp-folder temp \
    --alignments-per-thread 1000000 \
    --nthreads ${cpu}

    gzip gene-metrics.csv

  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD" 
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File gene_metrics = "gene-metrics.csv.gz"
  }
}
