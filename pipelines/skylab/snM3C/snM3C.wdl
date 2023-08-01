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
  # version of the pipeline
  String pipeline_version = "1.0.0"
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
      genome_fa = genome_fa,
      plate_id = plate_id
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

    String docker_image = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
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

    python3 <<CODE
    import re
    import os

    # Parsing stats.txt file
    stats_file_path = '/cromwell_root/~{plate_id}.stats.txt'
    adapter_counts = {}
    with open(stats_file_path, 'r') as file:
        content = file.read()

    adapter_matches = re.findall(r'=== First read: Adapter (\w+) ===\n\nSequence: .+; Type: .+; Length: \d+; Trimmed: (\d+) times', content)
    for adapter_match in adapter_matches:
        adapter_name = adapter_match[0]
        trimmed_count = int(adapter_match[1])
        adapter_counts[adapter_name] = trimmed_count

    # Removing fastq files with trimmed reads greater than 30
    directory_path = '/cromwell_root'
    threshold = 10000000

    for filename in os.listdir(directory_path):
        if filename.endswith('.fq.gz'):
            file_path = os.path.join(directory_path, filename)
            adapter_name = re.search(r'A(\d+)-R', filename)
            if adapter_name:
                adapter_name = 'A' + adapter_name.group(1)
                if adapter_name in adapter_counts and adapter_counts[adapter_name] > threshold:
                    os.remove(file_path)
                    print(f'Removed file: {filename}')
    CODE

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

    String docker_image = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
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
