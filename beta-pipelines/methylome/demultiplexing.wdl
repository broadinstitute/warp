version 1.0

workflow Demultiplexing {

  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_name
  }


  scatter(idx in range(length(fastq_input_read1))) {
   call Demultiplex {
    input:
      read1 = [fastq_input_read1[idx]],
      read2 = [fastq_input_read2[idx]],
      random_primer_indexes = random_primer_indexes,
      plate_name = plate_name
   }
  }
  output {
    Array[File] output_fastqs = flatten(Demultiplex.output_fastqs)
  }
}

  task Demultiplex {
    input {
      Array[File] read1
      Array[File] read2
      File random_primer_indexes
      String plate_name

      String docker_image = "ekiernan/yap_hisat:v4"
      Int disk_size = 50
      Int mem_size = 10
    }

  command <<<
    set -euo pipefail
    declare -a fastq1_file=(~{sep=' ' read1})
    declare -a fastq2_file=(~{sep=' ' read2})

    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels -g file:~{random_primer_indexes} -o ~{plate_name}.r1.{name}.output.fq.gz -p ~{plate_name}.r2.{name}.output.fq.gz $fastq1_file $fastq2_file > ~{plate_name}.stats.txt
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }


  output {
    Array[File] output_fastqs = glob("*.fq.gz")
    File stats = "~{plate_name}.stats.txt"
  }
 }
