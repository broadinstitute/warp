version 1.0

task CalculateGeneMetrics {
  input {
    File tsv_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:sctools-optimized"
    Int machine_mem_mb = 30000
    Int cpu = 1
    Int disk = ceil(size(tsv_input, "Gi") * 4)
    Int preemptible = 3
  }

  meta {
    description: "Calculate gene metrics from data in tsv_input."
  }

  parameter_meta {
    tsv_input: "a tsv files with sorted tags according to the order of GE, CB and UB."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    CalculateGeneMetricsFast -i "${tsv_input}" -o gene-metrics.csv.gz
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
    File tsv_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:sctools-optimized"
    Int machine_mem_mb = 45000
    Int cpu = 1
    Int disk = ceil(size(tsv_input, "Gi") * 2)
    Int preemptible = 3
  }

  meta {
    description: "Calculate cell metrics from data in tsv_input."
  }

  parameter_meta {
    tsv_input: "A tsv file augmented with CB, UB and GE tags."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    CalculateCellMetricsFast -i "~{tsv_input}" -a t -o cell-metrics.csv.gz
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
