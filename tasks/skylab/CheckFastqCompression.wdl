version 1.0

task CheckCompression {
  input {
    File fastq

    String FQ = "defaultFQ.gz"

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-hisat2:v0.4.0-fk"
  Int machine_mem_mb = 16500
  Int cpu = 4
  # Using (fastq) x 100 gives factor of a few buffer. BAM can be up to ~5 x (fastq1 + fastq2).
  # Need room for unsorted + sorted bam + temp sorting space + zipped and unzipped ref. Add 10 GiB buffer.
  Int disk = ceil((size(fastq, "GiB")) * 100 + 10)
  Int preemptible = 5
}
  meta {
    description: "HISAT2 alignment task will align paired-end fastq reads to reference genome."
  }

  parameter_meta {
    fastq: "gz forward fastq file"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    # Note that files MUST be gzipped or the module will not function properly
    # This will be addressed in the future either by a change in how Hisat2 functions or a more
    # robust test for compression type.

    # fix names if necessary.
    if (file "${fastq}" | grep -q compressed); then
        if [[ "${fastq}" != *.gz ]]; then
            if [[ "${fastq}" != *.fastq ]]; then
                FQ=${fastq}.fastq.gz
                mv ${fastq} ${fastq}.fastq.gz
            else
                FQ=${fastq}.gz
                mv ${fastq} ${fastq}.gz
            fi
        else
            FQ=${fastq}
        fi
    fi

  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    String fastq_name = "${FQ}"
  }
}
