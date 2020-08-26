version 1.0

task FastqProcessing {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    File whitelist

    # runtime values
    String docker = "quay.io/humancellatlas/fastq-process"
    Int machine_mem_mb = 3850
    Int cpu = 16   
    #TODO decided cpu
    # estimate that bam is approximately equal in size to fastq, add 20% buffer
    Int disk = 1200
    #TODO fix ceil(size(fastq_file, "GiB") * 2.2)

    Int preemptible = 3
  }

  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500

  meta {
    description: "Converts a fastq file into an unaligned bam file."
  }

  parameter_meta {
    r1_fastq: "input fastq file"
    r2_fastq: "input fastq file"
    i1_fastq: "input fastq file"
    whitelist: "10x genomics cell barcode whitelist"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    fastqprocess \
       --bam-size 1.0 \
       --I1 ${sep=' --I1 ' i1_fastq} \
       --R1 ${sep=' --R1 ' r1_fastq} \
       --R2 ${sep=' --R2 ' r2_fastq} \
       --white-list "${whitelist}"
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    Array[File] bam_output_array = glob("subfile_*")
  }
}
