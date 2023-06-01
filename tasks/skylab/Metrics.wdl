version 1.0

task CalculateCellMetrics {
  input {
    File bam_input
    File original_gtf
    File? mt_genes
    String input_id

    # runtime values

    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1679941323"
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
    --metric-output "~{input_id}.cell-metrics.csv" \
    --compute-metric \
    --metric-type cell \
    --barcode-tag CB \
    --umi-tag UB \
    --gene-tag GX \
    --temp-folder temp \
    --alignments-per-thread 1000000 \
    --nthreads ${cpu} \
    ~{"--mitochondrial-gene-names-filename " + mt_genes}

    # add the column header "CellID" to the first column in the .csv file
    sed '1s/^/CellID/' ~{input_id}.cell-metrics.csv > updated.~{input_id}.cell-metrics.csv

    # remove the following columns: reads_unmapped, reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, duplicate_reads, reads_mapped_intergenic
    cut -d',' -f 1-4,8-9,11-26,29-36 updated.~{input_id}.cell-metrics.csv > ~{input_id}.cell-metrics.csv

    gzip ~{input_id}.cell-metrics.csv
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
    File cell_metrics = "~{input_id}.cell-metrics.csv.gz"
  }
}

task CalculateGeneMetrics {
  input {
    File bam_input
    File? mt_genes
    String input_id
    # runtime values

    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1679941323"
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
    --metric-output "~{input_id}.gene-metrics.csv" \
    --compute-metric \
    --metric-type gene \
    --gene-tag GX \
    --barcode-tag CB \
    --umi-tag UB \
    --temp-folder temp \
    --alignments-per-thread 1000000 \
    --nthreads ${cpu} \
    ~{"--mitochondrial-gene-names-filename " + mt_genes}

    # add the column header "ID" to the first column in the .csv file
    sed '1s/^/GeneID/' ~{input_id}.gene-metrics.csv > updated.~{input_id}.gene-metrics.csv

    # remove the following columns: reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, duplicate_reads
    cut -d',' -f 1-4,8-9,11-27 updated.~{input_id}.gene-metrics.csv > ~{input_id}.gene-metrics.csv

    gzip ~{input_id}.gene-metrics.csv

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
    File gene_metrics = "~{input_id}.gene-metrics.csv.gz"
  }
}

task CalculateUMIsMetrics {
  input {
    File bam_input
    File? mt_genes
    String input_id
    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1679490798"
    Int machine_mem_mb = 16000
    Int cpu = 8
    Int disk = ceil(size(bam_input, "Gi") * 4)
    Int preemptible = 3
  }


  meta {
    description: "Sort bam_input by gene, then cell, then molecule."
  }

  parameter_meta {
    bam_input: "Input bam file containing reads marked with tags for cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
    input_id: "Name of sample matching this file; inserted into read group header"
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
    --metric-output "~{input_id}.umi-metrics.csv" \
    --compute-metric \
    --metric-type umi \
    --gene-tag GX \
    --barcode-tag CB \
    --umi-tag UB \
    --temp-folder temp \
    --alignments-per-thread 1000000 \
    --nthreads ${cpu} \
    ~{"--mitochondrial-gene-names-filename " + mt_genes}

    gzip ~{input_id}.umi-metrics.csv

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
    File umi_metrics = "~{input_id}.umi-metrics.csv.gz"
  }
}

task FastqMetricsSlideSeq {

  input {
    Array[File] r1_fastq
    String read_structure
    String sample_id
    File whitelist


    # Runtime attributes
    String docker =  "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1679490798"
    Int cpu = 16
    Int machine_mb = 40000
    Int disk = ceil(size(r1_fastq, "GiB")*3)  + 50
    Int preemptible = 3
  }

  meta {
    description: "Converts a set of fastq files to unaligned bam file, also corrects barcodes and partitions the alignments by barcodes. Allows for variable barcode and umi lengths as input"
  }

  parameter_meta {
    r1_fastq: "Array of Read 1 FASTQ files - forward read, contains cell barcodes and molecule barcodes"
    sample_id: "Name of sample matching this file, inserted into read group header"
    read_structure: "A string that specifies UMI (M) and Barcode (C) positions in the Read 1 fastq"
  }

  command <<<
    set -e

    cut -f 1 ~{whitelist} > WhiteList.txt
    declare -a arr_fastqs=(~{sep=' ' r1_fastq})
    p=" --R1 "
    arr_fastqs=( "${arr_fastqs[@]/#/$p}" ) 
    fastq_metrics \
    --white-list WhiteList.txt \
    --read-structure "~{read_structure}" \
    --sample-id "~{sample_id}" \
    ${arr_fastqs[@]}
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mb} MiB"
    disks: "local-disk ${disk} HDD"
    preemptible: preemptible
  }

  output {
    File barcode_distribution = "~{sample_id}.barcode_distribution_XC.txt"
    File umi_distribution = "~{sample_id}.barcode_distribution_XM.txt"
    File numReads_perCell = "~{sample_id}.numReads_perCell_XC.txt"
    File numReads_perUMI = "~{sample_id}.numReads_perCell_XM.txt"
  }
}

