version 1.0
import "../../../tasks/broad/Utilities.wdl" as utils


workflow snm3C {

    input {
        Array[File] fastq_input_read1
        Array[File] fastq_input_read2
        File random_primer_indexes
        String plate_id
        String cloud_provider
        # mapping inputs
        File tarred_index_files
        File genome_fa
        File chromosome_sizes

        String r1_adapter = "AGATCGGAAGAGCACACGTCTGAAC"
        String r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGA"
        Int r1_left_cut = 10
        Int r1_right_cut = 10
        Int r2_left_cut = 10
        Int r2_right_cut = 10
        Int min_read_length = 30
        Int num_upstr_bases = 0
        Int num_downstr_bases = 2
        Int compress_level = 5
        Int batch_number
    }
    #docker images
    String m3c_yap_hisat_docker = "m3c-yap-hisat:2.4"
    # Determine docker prefix based on cloud provider
    String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
    String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix
    String cromwell_root_dir = if cloud_provider == "gcp" then "/cromwell_root" else "/cromwell-executions"

    # make sure either gcp or azr is supplied as cloud_provider input
    if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
        call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
        input:
            message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
        }
    }

    # version of the pipeline
    String pipeline_version = "4.0.1"

    call Demultiplexing {
        input:
            fastq_input_read1 = fastq_input_read1,
            fastq_input_read2 = fastq_input_read2,
            random_primer_indexes = random_primer_indexes,
            plate_id = plate_id,
            batch_number = batch_number,
            docker = docker_prefix + m3c_yap_hisat_docker,
            cromwell_root_dir = cromwell_root_dir
    }

    output {
        File testOutput = Demultiplexing.testOutput
    }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    Int batch_number
    String docker
    String cromwell_root_dir

    Int disk_size = 1000
    Int mem_size = 10
    Int preemptible_tries = 2
    Int cpu = 8
  }

  command <<<
    echo "TEST"
    set -euo pipefail

    ls -lR
    pwd


    # Cat files for each r1, r2
    cat ~{sep=' ' fastq_input_read1} > ~{cromwell_root_dir}/r1.fastq.gz
    cat ~{sep=' ' fastq_input_read2} > ~{cromwell_root_dir}/r2.fastq.gz

    ls -lRh

    # Run cutadapt
    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels -j 8 \
    -g file:~{random_primer_indexes} \
    -o ~{plate_id}-{name}-R1.fq.gz \
    -p ~{plate_id}-{name}-R2.fq.gz \
    r1.fastq.gz \
    r2.fastq.gz \
    > ~{cromwell_root_dir}/~{plate_id}.stats.txt



    touch test.txt

  >>>

  runtime {
    docker: docker
    disks: "local-disk ${disk_size} SSD"
    cpu: cpu
    memory: "${mem_size} GiB"
    preemptible: preemptible_tries
  }

  output {
    File testOutput = "test.txt"
    File r2 =  "~{cromwell_root_dir}/r2.fastq.gz"
    File r1 = "~{cromwell_root_dir}/r1.fastq.gz"
    }
}
