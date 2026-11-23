version 1.0

task SplitBamByCellBarcode {
  input {
    Array[File] bams_to_split
    Float size_in_mb = 1024.0

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"

    Int machine_mem_mb = 15258
    Int cpu = 16

    # we can calculate disk size for arrays of input files in WDL 1.0
    Int disk = 500
    # by default request non preemptible machine to make sure the slow cell barcode split step completes
    Int preemptible = 0
  }

  meta {
    description: "Splits a bam file into chunks of size_in_mb, guaranteeing that all information for each cell is fully contained in only one of the chunks"
  }

  parameter_meta {
    bams_to_split: "input bam files to split by barcode"
    size_in_mb: "target size for each output chunk"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    SplitBam \
      --bamfile ${sep=' ' bams_to_split} \
      --output-prefix subfile \
      --subfile-size ${size_in_mb} \
      --tags CB CR \
      --num-processes ${cpu}
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
    Array[File] bam_output_array = glob("subfile_*")
  }
}
