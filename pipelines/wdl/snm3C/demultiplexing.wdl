version 1.0
import "../../../tasks/wdl/Utilities.wdl" as utils

workflow snm3C_demultiplexing_only {

    input {
        Array[File] fastq_input_read1
        Array[File] fastq_input_read2
        File random_primer_indexes
        String plate_id
        String cloud_provider
        # mapping inputs
        Int batch_number = 6
    }
    #docker images
    String m3c_yap_hisat_docker = "m3c-yap-hisat:2.4"
    # Determine docker prefix based on cloud provider
    String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
    String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix
    String cromwell_root_dir = if cloud_provider == "gcp" then "/mnt/disks/cromwell_root" else "/cromwell-executions"

    # make sure either gcp or azr is supplied as cloud_provider input
    if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
        call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
        input:
            message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
        }
    }

    # version of the pipeline
    String pipeline_version = "0.0.1"

    call Demultiplexing {
        input:
            fastq_input_read1 = fastq_input_read1,
            fastq_input_read2 = fastq_input_read2,
            random_primer_indexes = random_primer_indexes,
            plate_id = plate_id,
            batch_number = batch_number,
            docker = docker_prefix + m3c_yap_hisat_docker
    }

    meta {
        allowNestedInputs: true
    }

    output {
        Array[File] tarred_demultiplexed_fastqs = Demultiplexing.tarred_demultiplexed_fastqs
        File stats = Demultiplexing.stats
    }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    Int batch_number
    Int min_threshold = 100
    Int max_threshold = 6000000
    String docker

    Int disk_size = 1000
    Int mem_size = 10
    Int preemptible_tries = 2
    Int cpu = 8
  }

  command <<<
    set -euo pipefail
    WORKING_DIR=`pwd`

    # Cat files for each r1, r2
    cat ~{sep=' ' fastq_input_read1} > $WORKING_DIR/r1.fastq.gz
    cat ~{sep=' ' fastq_input_read2} > $WORKING_DIR/r2.fastq.gz

    # Run cutadapt
    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels -j 8 \
    -g file:~{random_primer_indexes} \
    -o ~{plate_id}-{name}-R1.fq.gz \
    -p ~{plate_id}-{name}-R2.fq.gz \
    $WORKING_DIR/r1.fastq.gz \
    $WORKING_DIR/r2.fastq.gz \
    > $WORKING_DIR/~{plate_id}.stats.txt

    # Remove the fastq files that end in unknown-R1.fq.gz and unknown-R2.fq.gz
    rm $WORKING_DIR/*-unknown-R{1,2}.fq.gz

    python3 <<CODE
    import re
    import os

    # Parsing stats.txt file
    working_dir = os.getcwd()
    stats_file_path = os.path.join(working_dir, '~{plate_id}.stats.txt')
    adapter_counts = {}
    with open(stats_file_path, 'r') as file:
        content = file.read()

    adapter_matches = re.findall(r'=== First read: Adapter (\w+) ===\n\nSequence: .+; Type: .+; Length: \d+; Trimmed: (\d+) times', content)
    for adapter_match in adapter_matches:
        adapter_name = adapter_match[0]
        trimmed_count = int(adapter_match[1])
        adapter_counts[adapter_name] = trimmed_count

    # Removing fastq files with trimmed reads greater than 10000000 or less than 100
    for filename in os.listdir(working_dir):
        if filename.endswith('.fq.gz'):
            file_path = os.path.join(working_dir, filename)
            adapter_name = re.search(r'([A-Za-z]\d+)-R', filename).group(1)
            if adapter_name:
                if adapter_name in adapter_counts:
                    if adapter_counts[adapter_name] < ~{min_threshold} or adapter_counts[adapter_name] > ~{max_threshold}:
                        print("Removing ", file_path, " with count equal to ", adapter_counts[adapter_name])
                        os.remove(file_path)
    CODE

    # Check if the number of *R1.fq.gz files is 0
    if [[ $(ls | grep "\-R1.fq.gz" | wc -l) -eq 0 ]]; then
        echo "Error: No files found. All fastq files were removed. Exiting."
        exit 1
    fi

    # Batch the fastq files into folders of batch_number size
    R1_files=($(ls $WORKING_DIR | grep "\-R1.fq.gz"))
    R2_files=($(ls $WORKING_DIR | grep "\-R2.fq.gz"))
    batch_number=~{batch_number}
    total_files=${#R1_files[@]}
    echo "Total files: $total_files"

    if [[ $total_files -lt $batch_number ]]; then
        echo "Warning: Number of files is less than the batch number. Updating batch number to $total_files."
        batch_number=$total_files
    fi

    for i in $(seq 1 "${batch_number}"); do  # Use seq for reliable brace expansion
        mkdir -p "batch${i}"  # Combine batch and i, use -p to create parent dirs
    done

    # Counter for the folder index and create emptycells file
    folder_index=1

    # Distribute the FASTQ files and create TAR files
    for file in "${R1_files[@]}"; do
        sample_id=$(basename "$file" "-R1.fq.gz")
        r2_file="${sample_id}-R2.fq.gz"

        mv $WORKING_DIR/$file batch$((folder_index))/$file
        mv $WORKING_DIR/$r2_file batch$((folder_index))/$r2_file
        # Increment the counter
        folder_index=$(( (folder_index % $batch_number) + 1 ))
    done

    # Tar up files per batch
    echo "TAR files"
    for i in $(seq 1 "${batch_number}"); do
        tar -cf - $WORKING_DIR/batch${i}/*.fq.gz | pigz > ~{plate_id}.${i}.cutadapt_output_files.tar.gz
    done
  >>>

  runtime {
    docker: docker
    disks: "local-disk ${disk_size} SSD"
    cpu: cpu
    memory: "${mem_size} GiB"
    preemptible: preemptible_tries
  }

  output {
    Array[File] tarred_demultiplexed_fastqs = glob("*.tar.gz")
    File stats = "~{plate_id}.stats.txt"
    }
}
