version 1.0

task FastqToUBam {
  input {
    File fastq_file
    String input_id

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    Int machine_mem_mb = 3850
    Int cpu = 1
    # estimate that bam is approximately equal in size to fastq, add 20% buffer
    Int disk = ceil(size(fastq_file, "GiB") * 2.2)
    Int preemptible = 3
  }

  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500

  meta {
    description: "Converts a fastq file into an unaligned bam file."
  }

  parameter_meta {
    fastq_file: "input fastq file"
    input_id: "name of sample matching this file, inserted into read group header"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {

    set -e

    if (file ~{fastq_file} | grep -q compressed); then
        if [[ ~{fastq_file} != *.gz ]]; then
            if [[ ~{fastq_file} != *.fastq ]]; then
                FQ="~{fastq_file}".fastq.gz
                mv  "~{fastq_file}" "~{fastq_file}".fastq.gz
            else
                FQ="~{fastq_file}".gz
                mv "~{fastq_file}" "~{fastq_file}".gz
            fi
        else
            FQ=~{fastq_file}
        fi
    elif [[ ~{fastq_file} != *.fastq ]]; then
      FQ="~{fastq_file}".fastq
      mv  "~{fastq_file}" "~{fastq_file}".fastq
    else
      FQ="~{fastq_file}"
    fi

    java -Xmx~{command_mem_mb}m -jar /usr/picard/picard.jar FastqToSam \
      FASTQ=$FQ \
      SORT_ORDER=unsorted \
      OUTPUT=bamfile.bam \
      SAMPLE_NAME="~{input_id}"
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
    File bam_output = "bamfile.bam"
  }
}
