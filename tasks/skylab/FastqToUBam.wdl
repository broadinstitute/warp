version 1.0

task FastqToUBam {
  input {
    File fastq_file
    String sample_id
    String fastq_suffix = ""

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
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
    sample_id: "name of sample matching this file, inserted into read group header"
    fastq_suffix: "a suffix to add to the fastq file; useful with mangled file IDs, since picard requires that the file end in .gz or it will not detect the gzipping."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    # Adds fastq_suffix if it is passed
    if [ ! -z "${fastq_suffix}" ];
    then
        mv "${fastq_file}" "${fastq_file}""${fastq_suffix}"
    fi

    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar FastqToSam \
      FASTQ="${fastq_file}""${fastq_suffix}" \
      SORT_ORDER=unsorted \
      OUTPUT=bamfile.bam \
      SAMPLE_NAME="${sample_id}"
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "bamfile.bam"
  }
}
