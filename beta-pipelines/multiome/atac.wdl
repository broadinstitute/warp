version 1.0

workflow ATAC {
  meta {
    description: "Processing for single-cell ATAC-seq data from the level of raw fastq reads. This is the first step of the multiome pipeline. ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used in molecular biology to assess genome-wide chromatin accessibility. This pipeline processes 10x Genomics Multiome ATAC FASTQ files."
    allowNestedInputs: true
  }

  input {
    # Fastq inputs
    File read1_fastq_gzipped
    File read2_fastq_gzipped
    File read3_fastq_gzipped

    # Output prefix/base name for all intermediate files and pipeline outputs
    String output_base_name
    
    # BWA ref 
    File tar_bwa_reference
    
    # script for monitoring tasks 
    File monitoring_script
  }

  parameter_meta {
    read1_fastq_gzipped: "read 1 FASTQ file as input for the pipeline, contains read 1 of paired reads"
    read2_fastq_gzipped: "read 2 FASTQ file as input for the pipeline, contains the cellular barcodes corresponding to the reads in the read1 FASTQ and read 3 FASTQ"
    read3_fastq_gzipped: "read 3 FASTQ file as input for the pipeline, contains read 2 of paired reads"
    output_base_name: "base name to be used for the pipelines output and intermediate files"
    monitoring_script : "script to monitor resource comsumption of tasks"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"

  }

  call AddBarcodes {
    input:
      read1_fastq = read1_fastq_gzipped,
      read3_fastq = read3_fastq_gzipped,
      barcodes_fastq = read2_fastq_gzipped,
      output_base_name = output_base_name
  }
  call TrimAdapters {
    input:
      fastq_input_read1 = AddBarcodes.fastq_barcodes_output_read1,
      fastq_input_read3 = AddBarcodes.fastq_barcodes_output_read3,
      output_base_name = output_base_name,
      monitoring_script = monitoring_script
   }
  call BWAPairedEndAlignment {
    input:
      fastq_input_read1 = TrimAdapters.fastq_trimmed_adapter_output_read1,
      fastq_input_read3 = TrimAdapters.fastq_trimmed_adapter_output_read3,
      tar_bwa_reference = tar_bwa_reference,
      output_base_name = output_base_name,
      monitoring_script = monitoring_script
  }
    
  output {
    File bam_aligned_output = BWAPairedEndAlignment.bam_aligned_output
  }   
}

  task AddBarcodes {
    input {
      File read1_fastq
      File read3_fastq
      File barcodes_fastq
      String output_base_name
      Int mem_size = 5
      String docker_image = "us.gcr.io/broad-gotc-prod/atac_barcodes:1.0.3-1679503564"
      Int disk_size = ceil(2 * ( size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(barcodes_fastq, "GiB") )) + 200
  }

   parameter_meta {
      read1_fastq: "Read 1 FASTQ with read 1 of paired reads"
      read3_fastq: "Read 3 FASTQ with read 2 of paired reads"
      barcodes_fastq: "Read 2 FASTQ with cellular barcodes"
      mem_size: "Size of memory in GB"
      output_base_name: "base name to be used for the output of the task"
      docker_image: "the docker image using cutadapt to be used (default: )"
      disk_size : "disk size used in trimming adapters step"
  }
      
    # output names for trimmed reads
    String fastq_barcodes_read1 = output_base_name + ".R1.barcodes.fastq"
    String fastq_barcodes_read3 = output_base_name + ".R3.barcodes.fastq"
   
    # Adding barcodes to read 1 and read 3 fastq
    command <<<
      set -euo pipefail
      mv ~{read1_fastq} r1.fastq.gz
      mv ~{read3_fastq} r3.fastq.gz
      mv ~{barcodes_fastq} barcodes.fastq.gz
      gunzip r1.fastq.gz r3.fastq.gz barcodes.fastq.gz
      python3 /usr/gitc/atac_barcodes.py -r1 r1.fastq -r3 r3.fastq -cb barcodes.fastq -out_r1 ~{fastq_barcodes_read1} -out_r3 ~{fastq_barcodes_read3}
      echo these are the zipped files and sizes
      ls -l
      gzip ~{fastq_barcodes_read1} ~{fastq_barcodes_read3}
      echo these are the zipped files
      ls -l
     >>>

    # use docker image for given tool cutadapat
    runtime {
      docker: docker_image
      disks: "local-disk ${disk_size} HDD"
      memory: "${mem_size} GiB"
    }

    output {
      File fastq_barcodes_output_read1 = "~{fastq_barcodes_read1}.gz"
      File fastq_barcodes_output_read3 = "~{fastq_barcodes_read3}.gz"
    }
  }
  # trim read 1 and read 2 adapter sequeunce with cutadapt
  task TrimAdapters {
    input {
      File fastq_input_read1
      File fastq_input_read3
      String output_base_name
      String docker_image = "quay.io/broadinstitute/cutadapt:1.18"
      File monitoring_script
      Int disk_size = ceil(2 * ( size(fastq_input_read1, "GiB") + size(fastq_input_read3, "GiB") )) + 200
      Int mem_size = 4
      Int min_length = 10
      Int quality_cutoff = 0
      String adapter_seq_read1
      String adapter_seq_read3
  }

   parameter_meta {
      fastq_input_read1: "read 1 fastq file containing sequencing reads as input for the pipeline"
      fastq_input_read3: "read 3 fastq file containing sequencing reads as input for the pipeline"
      min_length: "the minimum length for trimming. Reads that are too short even before adapter removal are also discarded"
      quality_cutoff: "cutadapt option to trim low-quality ends from reads before adapter removal"
      adapter_seq_read1: "cutadapt option for the sequence adapter for read 1 fastq"
      adapter_seq_read3: "cutadapt option for the sequence adapter for read 3 fastq"
      output_base_name: "base name to be used for the output of the task"
      docker_image: "the docker image using cutadapt to be used (default: quay.io/broadinstitute/cutadapt:1.18)"
      monitoring_script : "script to monitor resource consumption of tasks"
      mem_size: "the size of memory used during trimming adapters"
      disk_size : "disk size used in trimming adapters step"
  }
      
    # output names for trimmed reads
    String fastq_trimmed_adapter_output_name_read1 = output_base_name + ".R1.trimmed_adapters.fastq.gz"
    String fastq_trimmed_adapter_output_name_read3 = output_base_name + ".R3.trimmed_adapters.fastq.gz"
   
    # using cutadapt to trim off sequence adapters
    command <<<
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
        -A ~{adapter_seq_read3} \
        --output ~{fastq_trimmed_adapter_output_name_read1} \
        --paired-output ~{fastq_trimmed_adapter_output_name_read3} \
        ~{fastq_input_read1} ~{fastq_input_read3}
  >>>

    # use docker image for given tool cutadapat
    runtime {
      docker: docker_image
      disks: "local-disk ${disk_size} HDD"
      memory: "${mem_size} GiB" 
  }

    output {
      File fastq_trimmed_adapter_output_read1 = fastq_trimmed_adapter_output_name_read1
      File fastq_trimmed_adapter_output_read3 = fastq_trimmed_adapter_output_name_read3
      File monitoring_log = "monitoring.log"
    }
  }

  # align the two trimmed fastq as piared end data using BWA
  task BWAPairedEndAlignment {
    input {
      File fastq_input_read1
      File fastq_input_read3
      File tar_bwa_reference
      String read_group_id = "RG1"
      String read_group_sample_name = "RGSN1"
      String output_base_name
      String docker_image = "us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091"
      File monitoring_script
      Int disk_size = ceil(3.25 * (size(fastq_input_read1, "GiB") + size(fastq_input_read3, "GiB") + size(tar_bwa_reference, "GiB"))) + 200 
      Int nthreads = 16
      Int mem_size = 8
   }

    parameter_meta {
      fastq_input_read1: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
      fastq_input_read3: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
      tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
      read_group_id: "the read group id to be added upon alignment"
      read_group_sample_name: "the read group sample to be added upon alignment"
      nthreads: "the number of threads to use during bwa alignment"
      mem_size: "the size of memory used during alignment"
      disk_size : "disk size used in bwa alignment step"
      output_base_name: "basename to be used for the output of the task"
      docker_image: "the docker image using BWA to be used (default: us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091)"
      monitoring_script : "script to monitor resource comsumption of tasks"
    }

    String bam_aligned_output_name = output_base_name + ".aligned.bam"

    # bwa and call samtools to convert sam to bam
    command <<<

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
        ~{fastq_input_read1} ~{fastq_input_read3} \
        | samtools view -bS - > ~{bam_aligned_output_name}    
     >>>

    runtime {
      docker: docker_image
      disks: "local-disk ${disk_size} HDD"
      cpu: nthreads
      memory: "${mem_size} GiB" 
    }


    output {
      File bam_aligned_output = bam_aligned_output_name
      File monitoring_log = "monitoring.log"
    }
  }
