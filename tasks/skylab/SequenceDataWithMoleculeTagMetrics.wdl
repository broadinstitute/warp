version 1.0

task CalculateGeneMetrics {
  input {
    File bam_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 22000
    Int cpu = 1
    Int disk = ceil(size(bam_input, "Gi") * 4)
    Int preemptible = 3
  }

  meta {
    description: "Calculate gene metrics from the reads in bam_input."
  }

  parameter_meta {
    bam_input: "an aligned bam file augmented with CB, UB, GE, CY, UY, and XF tags."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    CalculateGeneMetrics -i "${bam_input}" -o gene-metrics.csv.gz
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

task CalculateCellMetrics {
  input {
    File bam_input
    File original_gtf

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 45000
    Int cpu = 1
    Int disk = ceil(size(bam_input, "Gi") * 2)
    Int preemptible = 3
  }

  meta {
    description: "Calculate cell metrics from the reads in bam_input."
  }

  parameter_meta {
    bam_input: "An aligned bam file augmented with CB, UB, GE, CY, UY, and XF tags."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    CalculateCellMetrics -i "${bam_input}" -o cell-metrics.csv.gz --gtf-annotation-file "${original_gtf}"
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


task MergeGeneMetrics {
  input {
    Array[File] metric_files

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 3850
    Int cpu = 1
    Int disk = 20
    Int preemptible = 3
  }

  meta {
    description: "Merge an array of gene metric files with the same metric categories and potentially overlapping sets of gene features"
  }

  parameter_meta {
    metric_files: "A set of metrics files, each measuring a potentially overlapping set of genes in disjoint sets of cells"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    MergeGeneMetrics -o merged-gene-metrics.csv.gz ${sep=' ' metric_files}
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File gene_metrics = "merged-gene-metrics.csv.gz"
  }
}

task MergeCellMetrics {
  input {
    Array[File] metric_files

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 3850
    Int cpu = 1
    Int disk = 20
    Int preemptible = 3
  }

  meta {
    description: "Concatenate multiple cell metrics files into a single matrix"
  }

  parameter_meta {
    metric_files: "An array of cell metrics files that contain the same metric types, but different sets of cells"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    MergeCellMetrics -o merged-cell-metrics.csv.gz ${sep=' ' metric_files}
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File cell_metrics = "merged-cell-metrics.csv.gz"
  }
}
