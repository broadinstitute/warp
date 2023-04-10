version 1.0

workflow Demultiplexing {

  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String well_id
  }


  scatter(idx in range(length(fastq_input_read1))) {
   call Demultiplex {
    input:
      fastq_input_read1 = [fastq_input_read1[idx]],
      fastq_input_read2 = [fastq_input_read2[idx]],
      random_primer_indexes = random_primer_indexes,
      well_id = well_id
   }
  }
  output {
    Array[File] output_fastqs = Demultiplex.output_fastqs
  }
}

  task Demultiplex {
    input {
      Array[File] fastq_input_read1
      Array[File] fastq_input_read2
      File random_primer_indexes
      String well_id

      String docker_image = "ekiernan/yap_hisat:v4"
      Int disk_size = 50
      Int mem_size = 10
    }

  command <<<
    set -euo pipefail
    declare -a fastq1_file=(~{sep=' ' fastq_input_read1})
    declare -a fastq2_file=(~{sep=' ' fastq_input_read2})

    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels \
    -g file:~{random_primer_indexes} \
    -o ~{well_id}-{name}-R1.fq.gz \
    -p ~{well_id}-{name}-R2.fq.gz \
    $fastq1_file \
    $fastq2_file \
    > ~{well_id}.stats.txt

    # remove the fastq files that end in unknown-R1.fq.gz and unknown-R2.fq.gz
    rm *-unknown-R{1,2}.fq.gz

    # zip up all the output fq.gz files
    tar -zcvf ~{well_id}.cutadapt_output_files.tar.gz *.fq.gz
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File output_fastqs = "~{well_id}.cutadapt_output_files.tar.gz"
    File stats = "~{well_id}.stats.txt"
  }
 }
