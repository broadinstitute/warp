version 1.0

task CellSortBam {
  input {
    File bam_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 100000
    Int cpu = 2
    Int disk = ceil(size(bam_input, "Gi") * 8)
    Int preemptible = 3
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

    TagSortBam -i ${bam_input} -o cell-sorted.bam -t CB UB GE
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "cell-sorted.bam"
  }
}

task GeneSortBam {
  input {
    File bam_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 100000
    Int cpu = 2
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

    TagSortBam -i ${bam_input} -o gene-sorted.bam -t GE CB UB
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "gene-sorted.bam"
  }
}
