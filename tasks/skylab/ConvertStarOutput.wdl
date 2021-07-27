version 1.0

task ConvertStarOutput {

  input {
    File barcodes
    File features
    File matrix

    #runtime values
    String docker = "quay.io/humancellatlas/snss2-featurecount:0.1.0"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(matrix, "Gi") * 2) + 10
    Int preemptible = 3
  }

  meta {
    description: "Counts the exonic and intronic reads from a bamfile using featureCounts."
  }

  parameter_meta {
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

   # create the  compresed raw count matrix with the counts, gene names and the barcodes
    python create-npz-output.py \
        --barcodes ~{barcodes} \
        --features ~{features} \
        --matrix ~{matrix}

  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File sparse_counts_row_index = "sparse_counts_row_index.npy"
    File sparse_counts_col_index = "sparse_counts_col_index.npy"
    File sparse_counts = "sparse_counts.npz"
  }
}
