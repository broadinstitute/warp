version 1.0

workflow ATAC {
  meta {
    description: "Processing for single-cell ATAC-seq data from the level of raw fastq reads. This is the first step of the multiome pipeline. ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used in molecular biology to assess genome-wide chromatin accessibility. This pipeline accepts fastq files where the cell barcode has been added to the fastq read names as the first field."
    allowNestedInputs: true
  }

  input {
    # Fastq inputs
    File fastq_gzipped_input_read1
    File fastq_gzipped_input_read2

    # Output prefix/base name for all intermediate files and pipeline outputs
    String output_base_name

    # script for monitoring tasks 
    File monitoring_script
  }

  parameter_meta {
    fastq_gzipped_input_read1: "read 1 fastq file as input for the pipeline, the cellular barcodes must be the first part of the read name seperated by colon"
    fastq_gzipped_input_read2: "read 2 fastq file as input for the pipeline, the cellular barcodes must be the first part of the read name separated by colon"
    output_base_name: "base name to be used for the pipelines output and intermediate files"
    monitoring_script : "script to monitor resource comsumption of tasks"
  }

  call TrimAdapters {
    input:
      fastq_input_read1 = fastq_gzipped_input_read1,
      fastq_input_read2 = fastq_gzipped_input_read2,
      output_base_name = output_base_name,
      monitoring_script = monitoring_script
    }

  call BWAPairedEndAlignment {
    input:
      fastq_input_read1 = TrimAdapters.fastq_trimmed_adapter_output_read1,
      fastq_input_read2 = TrimAdapters.fastq_trimmed_adapter_output_read2,
      output_base_name = output_base_name,
      monitoring_script = monitoring_script
    }
}

  # trim read 1 and read 2 adapter sequeunce with cutadapt
  task TrimAdapters {
    input {
      File fastq_input_read1
      File fastq_input_read2
      String output_base_name
      String docker_image = "quay.io/broadinstitute/cutadapt:1.18"
      File monitoring_script
      Int disk_size = ceil(2 * ( size(fastq_input_read1, "GiB") + size(fastq_input_read2, "GiB") )) + 200
      Int mem_size = 4
      Int min_length 
      Int quality_cutoff
      String adapter_seq_read1
      String adapter_seq_read2
  }

   parameter_meta {
      fastq_input_read1: "read 1 fastq file as input for the pipeline"
      fastq_input_read2: "read 2 fastq file as input for the pipeline"
      min_length: "the minimum legnth for trimming. Reads that are too short even before adapter removal are also discarded"
      quality_cutoff: "cutadapt option to trim low-quality ends from reads before adapter removal"
      adapter_seq_read1: "cutadapt option for the sequence adapter for read 1 fastq"
      adapter_seq_read2: "cutadapt option for the sequence adapter for read 2 fastq"
      output_base_name: "base name to be used for the output of the task"
      docker_image: "the docker image using cutadapt to be used (default: quay.io/broadinstitute/cutadapt:1.18)"
      monitoring_script : "script to monitor resource comsumption of tasks"
      mem_size: "the size of memory used during trimming adapters"
      disk_size : "disk size used in trimming adapters step"
  }
      
    # output names for trimmed reads
    String fastq_trimmed_adapter_output_name_read1 = output_base_name + ".R1.trimmed_adapters.fastq.gz"
    String fastq_trimmed_adapter_output_name_read2 = output_base_name + ".R2.trimmed_adapters.fastq.gz"
   
    # using cutadapt to trim off sequence adapters
    command {
      set -euo pipefail

      if [ ! -z "~{monitoring_script}" ]; then
        chmod a+x ~{monitoring_script}
        ~{monitoring_script} > monitoring.log &
      else
        echo "No monitoring script given as input" > monitoring.log &
      fi

      # fastq's, "-f", -A for paired adapters read 2"
      cutadapt \
        -f fastq \
        --minimum-length ~{min_length} \
        --quality-cutoff ~{quality_cutoff} \
        --adapter ~{adapter_seq_read1} \
        -A ~{adapter_seq_read2} \
        --output ~{fastq_trimmed_adapter_output_name_read1} \
        --paired-output ~{fastq_trimmed_adapter_output_name_read2} \
        ~{fastq_input_read1} ~{fastq_input_read2}
  }

    # use docker image for given tool cutadapat
    runtime {
      docker: docker_image
      disks: "local-disk ${disk_size} HDD"
      memory: "${mem_size} GiB" 
  }

    output {
      File fastq_trimmed_adapter_output_read1 = fastq_trimmed_adapter_output_name_read1
      File fastq_trimmed_adapter_output_read2 = fastq_trimmed_adapter_output_name_read2
      File monitoring_log = "monitoring.log"
    }
  }

  # align the two trimmed fastq as piared end data using BWA
  task BWAPairedEndAlignment {
    input {
      File fastq_input_read1
      File fastq_input_read2
      File tar_bwa_reference
      String read_group_id = "RG1"
      String read_group_sample_name = "RGSN1"
      String output_base_name
      String docker_image = "us.gcr.io/broad-gotc-prod/bwa:1.0.0-0.7.17-1660770463"
      File monitoring_script
      Int disk_size = ceil(3.25 * (size(fastq_input_read1, "GiB") + size(fastq_input_read2, "GiB") + size(tar_bwa_reference, "GiB"))) + 200 
      Int nthreads = 16
      Int mem_size = 8
   }

    parameter_meta {
      fastq_input_read1: "the trimmed read 1 fastq file as input for the aligner"
      fastq_input_read2: "the trimmed read 1 fastq file as input for the aligner"
      tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
      read_group_id: "the read group id to be added upon alignment"
      read_group_sample_name: "the read group sample to be added upon alignment"
      nthreads: "the number of threads to use during bwa alignment"
      mem_size: "the size of memory used during alignment"
      disk_size : "disk size used in bwa alignment step"
      output_base_name: "basename to be used for the output of the task"
      docker_image: "the docker image using BWA to be used (default: us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730)"
      monitoring_script : "script to monitor resource comsumption of tasks"
    }

    String sam_aligned_output_name = output_base_name + ".aligned.sam"

    # sort with samtools
    command {

      set -euo pipefail

      if [ ! -z "~{monitoring_script}" ]; then
        chmod a+x ~{monitoring_script}
         ~{monitoring_script} > monitoring.log &
      else
        echo "No monitoring script given as input" > monitoring.log &
      fi

      # prepare reference
      declare -r REF_DIR=$(mktemp -d genome_referenceXXXXXX)
      tar -xf "~{tar_bwa_reference}" -C $REF_DIR --strip-components 1
      rm "~{tar_bwa_reference}"

      # align w/ BWA: -t for number of cores
      bwa \
        mem \
        -R "@RG\tID:~{read_group_id}\tSM:~{read_group_sample_name}" \
        -t ~{nthreads} \
        $REF_DIR/genome.fa \
        ~{fastq_input_read1} ~{fastq_input_read2} \
        > ~{sam_aligned_output_name}
    }

    runtime {
      docker: docker_image
      disks: "local-disk ${disk_size} HDD"
      cpu: nthreads
      memory: "${mem_size} GiB" 
    }

    output {
      File sam_aligned_output = sam_aligned_output_name
      File monitoring_log = "monitoring.log"
    }
  }
