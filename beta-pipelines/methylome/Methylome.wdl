version 1.0

workflow Methylome {

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


  scatter(idx in range(length(fastq_input_read1))) {
    call Demultiplexing {
      input:
        fastq_input_read1 = [fastq_input_read1[idx]],
        fastq_input_read2 = [fastq_input_read2[idx]],
        random_primer_indexes = random_primer_indexes,
        plate_id = plate_id,
        output_basename = output_basename + "_shard" + idx
    }
  }

  call Mapping {
    input:
      tarred_demultiplexed_fastqs = Demultiplexing.tarred_demultiplexed_fastqs,
      tarred_index_files = tarred_index_files,
      mapping_yaml = mapping_yaml,
      snakefile = snakefile,
      chromosome_sizes = chromosome_sizes,
      genome_fa = genome_fa,
  }
  output {
    File MappingSummary = Mapping.mappingSummary
    Array[File] allcFiles = Mapping.allcFiles
    Array[File] allc_CGNFiles = Mapping.allc_CGNFiles
    Array[File] bamFiles = Mapping.bamFiles
    Array[File] detail_statsFiles = Mapping.detail_statsFiles
    Array[File] hicFiles = Mapping.hicFiles

  }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    String output_basename

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
    -o ~{output_basename}-{name}-R1.fq.gz \
    -p ~{output_basename}-{name}-R2.fq.gz \
    $fastq1_file \
    $fastq2_file \
    > ~{output_basename}.stats.txt

    # remove the fastq files that end in unknown-R1.fq.gz and unknown-R2.fq.gz
    rm *-unknown-R{1,2}.fq.gz

    # zip up all the output fq.gz files
    tar -zcvf ~{output_basename}.cutadapt_output_files.tar.gz *.fq.gz
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File tarred_demultiplexed_fastqs = "~{output_basename}.cutadapt_output_files.tar.gz"
    File stats = "~{output_basename}.stats.txt"
    Array[File] demultiplexed_fastq_files = glob("*.fq.gz")
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

    echo "here are all the files extracted from the tar files:"
    ls -l


    #THIS TAR COMMAND IS NOT WORKING
    #tar -zxvf ~{sep=' ' tarred_demultiplexed_fastqs}


    # run the snakemake command
    cd ../
    echo "The snakemake command is being run here:"
    pwd

    /opt/conda/bin/snakemake --configfile mapping.yaml -j

    #ls a bunch of things so we can figure out what to grab as output of this task
    ls -l
    echo "ls /cromwell_root/group0/allc"
    ls -l /cromwell_root/group0/allc
    echo "ls /cromwell_root/group0/allc-CGN"
    ls -l /cromwell_root/group0/allc-CGN
    echo "ls /cromwell_root/group0/bam"
    ls -l /cromwell_root/group0/bam
    echo "ls /cromwell_root/group0/detail_stats"
    ls -l /cromwell_root/group0/detail_stats
    echo "ls /cromwell_root/group0/hic"
    ls -l /cromwell_root/group0/hic

  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File mappingSummary = "/cromwell_root/group0/MappingSummary.csv.gz"
    Array[File] allcFiles = glob("/cromwell_root/group0/allc/*")
    Array[File] allc_CGNFiles = glob("/cromwell_root/group0/allc-CGN/*")
    Array[File] bamFiles = glob("/cromwell_root/group0/bam/*")
    Array[File] detail_statsFiles = glob("/cromwell_root/group0/detail_stats/*")
    Array[File] hicFiles = glob("/cromwell_root/group0/hic/*")

  }
}
