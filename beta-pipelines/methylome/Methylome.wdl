version 1.0

workflow Methylome {

  input {
    # demultiplexing inputs
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id

    # mapping inputs
    File tarred_index_files
    File mapping_yaml
    File snakefile
    File chromosome_sizes
    File genome_fa
  }


  scatter(idx in range(length(fastq_input_read1))) {
    call Demultiplexing {
      input:
        fastq_input_read1 = [fastq_input_read1[idx]],
        fastq_input_read2 = [fastq_input_read2[idx]],
        random_primer_indexes = random_primer_indexes,
        plate_id = plate_id
    }
  }

  call Mapping {
    input:
      tarred_demultiplexed_fastqs = Demultiplexing.output_fastqs,
      tarred_index_files = tarred_index_files,
      mapping_yaml = mapping_yaml,
      snakefile = snakefile,
      chromosome_sizes = chromosome_sizes,
      genome_fa = genome_fa
    }

  output {
    Array[File] output_fastqs = Demultiplexing.output_fastqs
  }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id

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
    -o ~{plate_id}-{name}-R1.fq.gz \
    -p ~{plate_id}-{name}-R2.fq.gz \
    $fastq1_file \
    $fastq2_file \
    > ~{plate_id}.stats.txt

    # remove the fastq files that end in unknown-R1.fq.gz and unknown-R2.fq.gz
    rm *-unknown-R{1,2}.fq.gz

    # zip up all the output fq.gz files
    tar -zcvf ~{plate_id}.cutadapt_output_files.tar.gz *.fq.gz
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File output_fastqs = "~{plate_id}.cutadapt_output_files.tar.gz"
    File stats = "~{plate_id}.stats.txt"
  }
}


task Mapping {
  input {
    Array[File] tarred_demultiplexed_fastqs
    File tarred_index_files
    File mapping_yaml
    File snakefile
    File chromosome_sizes
    File genome_fa


    String docker_image = "ekiernan/yap_hisat:v4"
    Int disk_size = 200
    Int mem_size = 500
  }

  command <<<
    set -euo pipefail

    echo "pwd is"
    pwd
    echo "ls is"
    ls

    mkdir /group0
    mkdir /group0/reference/
    mkdir /group0/fastq/


    cp ~{tarred_index_files} /group0/reference/
    cp ~{chromosome_sizes} /group0/reference/
    cp ~{genome_fa} /group0/reference/
    cp ~{sep=' ' tarred_demultiplexed_fastqs} /group0/fastq/
    cp ~{mapping_yaml} /group0/
    cp ~{snakefile} /group0/



    # untar the index files
    cd /group0/reference/
    echo "Untarring the index files"
    tar -zxvf ~{tarred_index_files}
    rm ~{tarred_index_files}
    echo "The current working directory is (for the reference dir):"
    pwd
    echo "here is the ls command (for the reference dir):"
    ls



    # untar the demultiplexed fastq files
    cd ../fastq/
    echo "Untarring the demultiplexed fastq files"
    tar -zxvf ~{sep=' ' tarred_demultiplexed_fastqs}


    # run the snakemake command
    cd ../
    echo "The current working directory is  (the snakemake command is being run here:"
    pwd
    echo "here is the ls command (for the snakemake command):"
    ls

    /opt/conda/bin/snakemake --configfile mapping.yaml -j

  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    Array[File] unsorted_bam_files = "output.bam"
  }
}
