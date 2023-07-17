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

    # script for monitoring tasks
    File monitoring_script
  }

  call Demultiplexing {
    input:
      fastq_input_read1 = fastq_input_read1,
      fastq_input_read2 = fastq_input_read2,
      random_primer_indexes = random_primer_indexes,
      plate_id = plate_id,
      monitoring_script = monitoring_script
    }

  call Mapping {
    input:
      tarred_demultiplexed_fastqs = Demultiplexing.tarred_demultiplexed_fastqs,
      tarred_index_files = tarred_index_files,
      mapping_yaml = mapping_yaml,
      snakefile = snakefile,
      chromosome_sizes = chromosome_sizes,
      genome_fa = genome_fa,
      plate_id = plate_id,
      monitoring_script = monitoring_script
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

    String docker_image = "nikellepetrillo/yap-hisat:v8"
    Int disk_size = 50
    Int mem_size = 10

    # script for monitoring tasks
    File monitoring_script
  }

  command <<<
    set -euo pipefail

    if [ ! -z "~{monitoring_script}" ]; then
    chmod a+x ~{monitoring_script}
    ~{monitoring_script} > monitoring.log &
    else
    echo "No monitoring script given as input" > monitoring.log &
    fi

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

    # Count the number of reads in each fastq file and remove if over 10,000,000 reads. Also, remove its mate.
    for file in ~{plate_id}-*.fq.gz; do
      if [ -f $file ]; then
        num_reads=$(($(cat $file | wc -l) / 4))
        if [ $num_reads -gt 30 ]; then
          echo "Removing $file with $num_reads reads"
          rm $file

         # Remove the mate fastq file
          mate_file=${file/-R1./-R2.}
          if [ -f $mate_file ]; then
            echo "Removing the first $mate_file"
            rm $mate_file
          else
            mate_file=${file/-R2./-R1.}
            if [ -f $mate_file ]; then
              echo "Removing $mate_file"
              rm $mate_file
            fi
          fi
        fi
      fi
    done

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
    String plate_id
    # script for monitoring tasks
    File monitoring_script

    String docker_image = "nikellepetrillo/yap-hisat:v8"
    Int disk_size = 200
    Int mem_size = 500
  }

  command <<<
    set -euo pipefail

    if [ ! -z "~{monitoring_script}" ]; then
    chmod a+x ~{monitoring_script}
    ~{monitoring_script} > monitoring.log &
    else
    echo "No monitoring script given as input" > monitoring.log &
    fi

    mkdir group0/
    mkdir group0/fastq/
    mkdir group0/reference/

    cp ~{tarred_index_files} group0/reference/
    cp ~{chromosome_sizes} group0/reference/
    cp ~{genome_fa} group0/reference/
    cp ~{tarred_demultiplexed_fastqs} group0/fastq/
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
    tar -zxvf ~{tarred_demultiplexed_fastqs}
    rm ~{tarred_demultiplexed_fastqs}

    # run the snakemake command
    cd ../
    /opt/conda/bin/snakemake --configfile mapping.yaml -j

    # move outputs into /cromwell_root/
    mv /cromwell_root/group0/MappingSummary.csv.gz /cromwell_root/~{plate_id}_MappingSummary.csv.gz

    cd /cromwell_root/group0/allc
    tar -zcvf ~{plate_id}_allc_files.tar.gz *
    mv ~{plate_id}_allc_files.tar.gz /cromwell_root/
    cd ../allc-CGN
    tar -zcvf ~{plate_id}_allc-CGN_files.tar.gz *
    mv ~{plate_id}_allc-CGN_files.tar.gz /cromwell_root/
    cd ../bam
    tar -zcvf ~{plate_id}_bam_files.tar.gz *
    mv ~{plate_id}_bam_files.tar.gz /cromwell_root/
    cd ../detail_stats
    tar -zcvf ~{plate_id}_detail_stats_files.tar.gz *
    mv ~{plate_id}_detail_stats_files.tar.gz /cromwell_root/
    cd ../hic
    tar -zcvf ~{plate_id}_hic_files.tar.gz *
    mv ~{plate_id}_hic_files.tar.gz /cromwell_root/

  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File mappingSummary = "~{plate_id}_MappingSummary.csv.gz"
    File allcFiles = "~{plate_id}_allc_files.tar.gz"
    File allc_CGNFiles = "~{plate_id}_allc-CGN_files.tar.gz"
    File bamFiles = "~{plate_id}_bam_files.tar.gz"
    File detail_statsFiles = "~{plate_id}_detail_stats_files.tar.gz"
    File hicFiles = "~{plate_id}_hic_files.tar.gz"
  }
}