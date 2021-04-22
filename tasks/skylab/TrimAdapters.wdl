version 1.0

task TrimAdapters {

  input {
    File fastq1
    File? fastq2
    File adapter_list

    #runtime values
    String docker = "quay.io/humancellatlas/snss2-trim-adapters:0.1.0"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(fastq1, "Gi") * 2) + 10
    Int preemptible = 3
  }

  meta {
    description: "Trims adapters from FASTQ files."
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

    fastq-mcf \
       -C 200000 ~{adapter_list} \
       ~{fastq1} \
       ~{fastq2} \
       -o fastq_R1.trimmed.fastq.gz \
       -o fastq_R2.trimmed.fastq.gz
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File trimmed_fastq1 = "fastq_R1.trimmed.fastq.gz"
    File trimmed_fastq2 = "fastq_R2.trimmed.fastq.gz"
  }
}
