version 1.0

task CheckCompression {
  input {
    File fastq1
    File? fastq2

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
  Int machine_mem_mb = 16500
  Int cpu = 4
  # Using (fastq1 + fastq2) x 100 gives factor of a few buffer. BAM can be up to ~5 x (fastq1 + fastq2).
  # Need room for unsorted + sorted bam + temp sorting space + zipped and unzipped ref. Add 10 GiB buffer.
  Int disk = ceil((size(fastq1, "GiB") + size(fastq2, "GiB")) * 100 + size(hisat2_ref, "GiB") * 2 + 10)
  Int preemptible = 5
}
  meta {
    description: "HISAT2 alignment task will align paired-end fastq reads to reference genome."
  }

  parameter_meta {
    fastq1: "gz forward fastq file"
    fastq2: "gz reverse fastq file"
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
    if (file "${fastq1}" | grep -q compressed); then
        if [[ "${fastq1}" != *.gz ]]; then
            if [[ "${fastq1}" != *.fastq ]]; then
                FQ1=${fastq1}.fastq.gz
                mv ${fastq1} ${fastq1}.fastq.gz
            else
                FQ1=${fastq1}.gz
                mv ${fastq1} ${fastq1}.gz
            fi
        else
            FQ1=${fastq1}
        fi
    fi

    if (file "${fastq2}" | grep -q compressed); then
        if [[ "${fastq2}" != *.gz ]]; then
            if [[ "${fastq2}" != *.fastq ]]; then
                FQ2=${fastq2}.fastq.gz
                mv ${fastq2} ${fastq2}.fastq.gz
            else
                FQ2=${fastq2}.gz
                mv ${fastq2} ${fastq2}.gz
            fi
        else
            FQ2=${fastq2}
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
    File fastq1 = "${FQ1}"
    File? fastq2 = "${FQ2}"
  }
}
