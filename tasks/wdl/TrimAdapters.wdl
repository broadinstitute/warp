version 1.0

task TrimAdapters {

  input {
    Array[File] fastq1_input_files
    Array[File] fastq2_input_files
    File adapter_list
    Array[String] input_ids

    #runtime values
    String docker = "us.gcr.io/broad-gotc-prod/ea-utils:1.0.0-1.04.807-1659990665"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(2*(size(fastq1_input_files, "Gi") + size(fastq2_input_files, "Gi"))) + 10
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

  command <<<
    set -e

    declare -a fastq1_files=(~{sep=' ' fastq1_input_files})
    declare -a fastq2_files=(~{sep=' ' fastq2_input_files})
    declare -a output_prefix=(~{sep=' ' input_ids})
    for ((i=0; i<${#fastq1_files[@]}; ++i));
      do
        fastq1=${fastq1_files[$i]}
        fastq2=${fastq2_files[$i]}

        fastq-mcf \
           -C 200000 ~{adapter_list} \
           $fastq1 \
           $fastq2 \
           -o "${output_prefix[$i]}.trimmed_R1.fastq.gz" \
           -o "${output_prefix[$i]}.trimmed_R2.fastq.gz"
      done;
  >>>

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    Array[File] trimmed_fastq1_files = glob("*trimmed_R1.fastq.gz")
    Array[File] trimmed_fastq2_files = glob("*trimmed_R2.fastq.gz")
  }
}
