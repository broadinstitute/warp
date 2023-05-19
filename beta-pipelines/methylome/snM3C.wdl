version 1.0

workflow snM3C {

  input {
    # demultiplexing inputs
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    String output_basename = plate_id

    # mapping inputs
    File tarred_index_files
    File mapping_yaml
    File snakefile
    File chromosome_sizes
    File genome_fa
  }


  call Demultiplexing {
    input:
      fastq_input_read1 = fastq_input_read1,
      fastq_input_read2 = fastq_input_read2,
      random_primer_indexes = random_primer_indexes,
      plate_id = plate_id
    }

  call Mapping {
    input:
      tarred_demultiplexed_fastqs = Demultiplexing.tarred_demultiplexed_fastqs,
      tarred_index_files = tarred_index_files,
      mapping_yaml = mapping_yaml,
      snakefile = snakefile,
      chromosome_sizes = chromosome_sizes,
      genome_fa = genome_fa
    }

  output {
    File MappingSummary = Mapping.mappingSummary
    File allcFiles = Mapping.allcFiles
    File allc_CGNFiles = Mapping.allc_CGNFiles
    File bamFiles = Mapping.bamFiles
    File detail_statsFiles = Mapping.detail_statsFiles
    File hicFiles = Mapping.hicFiles
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

    # Cat files for each r1, r2
    cat ~{sep=' ' fastq_input_read1} > r1.fastq.gz
    cat ~{sep=' ' fastq_input_read2} > r2.fastq.gz

    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels \
    -g file:~{random_primer_indexes} \
    -o ~{plate_id}-{name}-R1.fq.gz \
    -p ~{plate_id}-{name}-R2.fq.gz \
    r1.fastq.gz \
    r2.fastq.gz \
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
    File tarred_demultiplexed_fastqs = "~{plate_id}.cutadapt_output_files.tar.gz"
    File stats = "~{plate_id}.stats.txt"}
}


task Mapping {
  input {
    File tarred_demultiplexed_fastqs
    File tarred_index_files
    File mapping_yaml
    File snakefile
    File chromosome_sizes
    File genome_fa

    String docker_image = "nikellepetrillo/yap-hisat:v8"
    Int disk_size = 200
    Int mem_size = 500
  }

  command <<<
    set -euo pipefail

    mkdir group0/
    mkdir group0/fastq/
    mkdir group0/reference/

    cp ~{tarred_index_files} group0/reference/
    cp ~{chromosome_sizes} group0/reference/
    cp ~{genome_fa} group0/reference/
    echo "copy tarred demulitplexed files to group0/fastq/"
    cp ~{sep=' ' tarred_demultiplexed_fastqs} group0/fastq/
    cp ~{mapping_yaml} group0/
    cp ~{snakefile} group0/

    # untar the index files
    cd group0/reference/
    echo "Untarring the index files"
    tar -zxvf ~{tarred_index_files}
    rm ~{tarred_index_files}
    samtools faidx hg38.fa

    # untar the demultiplexed fastq files
    cd ../fastq/

    # create array of .tar.gz files
    tar_files=(~{sep=' ' tarred_demultiplexed_fastqs})

    # loop through the array and untar all the tar files
    for file in "${tar_files[@]}"; do
        # Untar the file
         echo "Untarring $file"
         tar -xzf "$file"
    done

    # remove the tar files
    rm *tar.gz

    # run the snakemake command
    cd ../
    echo "The snakemake command is being run here:"
    pwd

    /opt/conda/bin/snakemake --configfile mapping.yaml -j

    mv /cromwell_root/group0/MappingSummary.csv.gz /cromwell_root/

    cd /cromwell_root/group0/allc
    tar -zcvf allc_files.tar.gz *
    mv allc_files.tar.gz /cromwell_root/
    cd ../allc-CGN
    tar -zcvf allc-CGN_files.tar.gz *
    mv allc-CGN_files.tar.gz /cromwell_root/
    cd ../bam
    tar -zcvf bam_files.tar.gz *
    mv bam_files.tar.gz /cromwell_root/
    cd ../detail_stats
    tar -zcvf detail_stats_files.tar.gz *
    mv detail_stats_files.tar.gz /cromwell_root/
    cd ../hic
    tar -zcvf hic_files.tar.gz *
    mv hic_files.tar.gz /cromwell_root/

  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File mappingSummary = "MappingSummary.csv.gz"
    File allcFiles = "allc_files.tar.gz"
    File allc_CGNFiles = "allc-CGN_files.tar.gz"
    File bamFiles = "bam_files.tar.gz"
    File detail_statsFiles = "detail_stats_files.tar.gz"
    File hicFiles = "hic_files.tar.gz"

  }
}
