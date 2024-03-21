version 1.0

task MergeSortBamFiles {
  input {
    Array[File] bam_inputs
    String sort_order
    String output_bam_filename

    Int compression_level = 5

    # runtime values
    String picard_cloud_docker_path
    Int machine_mem_mb = 18150
    Int cpu = 1
    # default to 500GiB of space
    Int disk = 500
    # by default request non preemptible machine to make sure the slow mergsort step completes
    Int preemptible = 0
  }

  # give the command 500MiB of overhead
  Int command_mem_mb = machine_mem_mb - 500

  meta {
    description: "Merge multiple bam files in the specified sort order"
  }

  parameter_meta {
    bam_inputs: "Merges Sam/Bam files"
    sort_order: "sort order of output bam"
    picard_cloud_docker_path: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    java -Dsamjdk.compression_level=${compression_level} -Xms${command_mem_mb}m -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar \
      MergeSamFiles \
      USE_THREADING=true \
      SORT_ORDER=${sort_order} \
      INPUT=${sep=' INPUT=' bam_inputs} \
      OUTPUT=~{output_bam_filename} \
  }

  runtime {
    docker: picard_cloud_docker_path
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }
  output {
    File output_bam = output_bam_filename
  }
}


