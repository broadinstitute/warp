version 1.0

workflow sampleFASTQ {
  meta {
    description: "Sample GEX data."
    allowNestedInputs: true
  }

  input {
    # Fastq inputs
    Array[File] read1_fastq_gzipped
    Array[File] read2_fastq_gzipped

    # Whitelist
    File whitelist
    
    # Read structure
    String read_structure

  }

  parameter_meta {
    read1_fastq_gzipped: "read 1 FASTQ file as input for the pipeline, contains read 1 of paired reads"
    read2_fastq_gzipped: "read 2 FASTQ file as input for the pipeline, contains the cellular barcodes corresponding to the reads in the read1 FASTQ and read 3 FASTQ"
    output_base_name: "base name to be used for the pipelines output and intermediate files"
  }

  scatter(idx in range(length(read1_fastq_gzipped))) {
    call SampleFastq {
      input:
        read1_fastq = read1_fastq_gzipped[idx],
        read2_fastq = read2_fastq_gzipped[idx],
        whitelist = whitelist,
        read_structure = read_structure
     }
  }
  
  output {
    Array[File] output_r1 = SampleFastq.fastq_R1_output
    Array[File] output_r2 = SampleFastq.fastq_R2_output
  }
}

task SampleFastq {

  input {
    File read1_fastq
    File read2_fastq
    String read_structure
    #String barcode_orientation = "FIRST_BP_RC"
    String output_base_name = basename(read1_fastq, ".fastq.gz")
    File whitelist

    # [?] copied from corresponding optimus wdl for fastqprocessing 
    # using the latest build of warp-tools in GCR
    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1690997141"
    # Runtime attributes [?]
    Int mem_size = 5
    Int cpu = 16   
    # TODO decided cpu
    # estimate that bam is approximately equal in size to fastq, add 20% buffer
    Int disk_size = ceil(2 * ( size(read1_fastq, "GiB") + size(read2_fastq, "GiB") )) + 400
    Int preemptible = 3
   }

  meta {
    description: "Converts a set of fastq files to unaligned bam file/fastq file, also corrects barcodes and partitions the alignments by barcodes. Allows for variable barcode and umi lengths as input, if applicable."
   }

  parameter_meta {
        read1_fastq: "Array of read 1 FASTQ files which contains the cellular barcodes"
		read2_fastq: "Array of read 2 FASTQ files of paired reads -- forward reads"
        output_base_name: "Name of sample matching this file, inserted into read group header"
        read_structure: "A string that specifies the barcode (C) and UMI positions in the Read 1 fastq"
        whitelist: "10x genomics cell barcode whitelist for GEX."
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
   }

  command <<<
  
    set -euo pipefail
    
    # Cat files for each r1 and r2 together
    cat ~{sep=' ' read1_fastq} > r1.fastq.gz
    cat ~{sep=' ' read2_fastq} > r2.fastq.gz
    
    # Call samplefastq
    # outputs fastq files with reads that have valid barcodes
    samplefastq \
        --sample-id "~{output_base_name}" \
        --R1 r1.fastq.gz \
        --R2 r2.fastq.gz \
        --white-list "~{whitelist}" \
        --output-format "FASTQ" \
        --read-structure "~{read_structure}"
        
        mv fastq_R1_0.fastq.gz "~{output_base_name}_fastq_R1_0.fastq.gz"
        mv fastq_R2_0.fastq.gz "~{output_base_name}_fastq_R2_0.fastq.gz"
        
    >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem_size} MiB"
    disks: "local-disk ${disk_size} HDD"
    preemptible: preemptible
   }
  
  output {
    File fastq_R1_output = "~{output_base_name}_fastq_R1_0.fastq.gz"
    File fastq_R2_output = "~{output_base_name}_fastq_R2_0.fastq.gz"
   }
 }



