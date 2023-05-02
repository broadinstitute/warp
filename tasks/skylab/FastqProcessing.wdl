version 1.0

task FastqProcessing {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    File whitelist
    Int chemistry
    String sample_id

    #using the latest build of warp-tools in GCR
    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1679490798"
    #runtime values
    Int machine_mem_mb = 40000
    Int cpu = 16   
    #TODO decided cpu
    # estimate that bam is approximately equal in size to fastq, add 20% buffer
    Int disk = ceil(size(r1_fastq, "GiB")*3 + size(r2_fastq, "GiB")*3) + 500

    Int preemptible = 3
  }

  meta {
    description: "Converts a set of fastq files to unaligned bam file, also corrects barcodes and partitions the alignments by barcodes."
  }

  parameter_meta {
    r1_fastq: "input fastq file"
    r2_fastq: "input fastq file"
    i1_fastq: "(optional) input fastq file"
    whitelist: "10x genomics cell barcode whitelist"
    chemistry: "chemistry employed, currently can be tenX_v2 or tenX_v3, the latter implies NO feature barcodes"
    sample_id: "name of sample matching this file, inserted into read group header"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    FASTQS=$(python3 <<CODE
    def rename_file(filename):
        import shutil
        import gzip
        import re
         
        iscompressed = True
        with gzip.open(filename, 'rt') as fin:
           try:
               _ = fin.readline()
           except:
               iscompressed = False

        basename = re.sub(r'.gz$', '', filename)
        basename = re.sub(r'.fastq$', '', basename)
 
        if iscompressed:
            # if it is already compressed then add an extension .fastq.gz
            newname = basename + ".fastq.gz" 
        else: 
            # otherwise, add just the .fastq extension
            newname = basename + ".fastq"

        if filename != newname:
            # safe to rename since the old and the new names are different
            shutil.move(filename, newname)

        return newname
    optstring = ""
     
    r1_fastqs = [ "${sep='", "' r1_fastq}" ]
    r2_fastqs = [ "${sep='", "' r2_fastq}" ]
    i1_fastqs = [ "${sep='", "' i1_fastq}" ]
    for fastq in r1_fastqs:
        if fastq.strip(): 
            optstring += " --R1 " + rename_file(fastq)
    for fastq in r2_fastqs:
        if fastq.strip(): 
            optstring += " --R2 " + rename_file(fastq)
    for fastq in i1_fastqs:
        if fastq.strip(): 
            optstring += " --I1 " + rename_file(fastq)
    print(optstring)
    CODE)

    # use the right UMI length depending on the chemistry
    if [ "~{chemistry}" == "2" ]; then
        ## V2
        UMILENGTH=10
    elif [ "~{chemistry}" == "3" ]; then
        ## V3
        UMILENGTH=12
    else
        echo Error: unknown chemistry value: "~{chemistry}"
        exit 1;
    fi

    fastqprocess \
        --bam-size 30.0 \
        --barcode-length 16 \
        --umi-length $UMILENGTH \
        --sample-id "~{sample_id}" \
        $FASTQS \
        --white-list "~{whitelist}" \
        --output-format FASTQ
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    Array[File] fastq_R1_output_array = glob("fastq_R1_*")
    Array[File] fastq_R2_output_array = glob("fastq_R2_*")
  }
}

task FastqProcessingSlidSeq {

  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    String read_structure
    String sample_id
    File whitelist


    # Runtime attributes
    String docker =  "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.14-test2"
    Int cpu = 16
    Int machine_mb = 40000
    Int disk = ceil(size(r1_fastq, "GiB")*3 + size(r2_fastq, "GiB")*3) + 50
    Int preemptible = 3
  }

  meta {
    description: "Converts a set of fastq files to unaligned bam file, also corrects barcodes and partitions the alignments by barcodes. Allows for variable barcode and umi lengths as input"
  }

  parameter_meta {
        r1_fastq: "Array of Read 1 FASTQ files - forward read, contains cell barcodes and molecule barcodes"
        r2_fastq: "Array of Read 2 FASTQ files - reverse read, contains cDNA fragment generated from captured mRNA"
        i1_fastq: "(optional) Array of i1 FASTQ files - index read, for demultiplexing of multiple samples on one flow cell."
        sample_id: "Name of sample matching this file, inserted into read group header"
        read_structure: "A string that specifies UMI (M) and Barcode (C) positions in the Read 1 fastq"

  }

  command {

    command {
    set -e

    FASTQS=$(python3 <<CODE
    def rename_file(filename):
        import shutil
        import gzip
        import re

        iscompressed = True
        with gzip.open(filename, 'rt') as fin:
          try:
              _ = fin.readline()
          except:
              iscompressed = False

        basename = re.sub(r'.gz$', '', filename)
        basename = re.sub(r'.fastq$', '', basename)

        if iscompressed:
            # if it is already compressed then add an extension .fastq.gz
            newname = basename + ".fastq.gz"
        else:
            # otherwis, add just the .fastq extension
            newname = basename + ".fastq"

        if filename != newname:
            # safe to rename since the old and the new names are different
            shutil.move(filename, newname)

        return newname
    optstring = ""

    r1_fastqs = [ "${sep='", "' r1_fastq}" ]
    r2_fastqs = [ "${sep='", "' r2_fastq}" ]
    i1_fastqs = [ "${sep='", "' i1_fastq}" ]
    for fastq in r1_fastqs:
        if fastq.strip():
            optstring += " --R1 " + rename_file(fastq)
    for fastq in r2_fastqs:
        if fastq.strip():
            optstring += " --R2 " + rename_file(fastq)
    for fastq in i1_fastqs:
        if fastq.strip():
            optstring += " --I1 " + rename_file(fastq)
    print(optstring)
    CODE)
    cut -f 1 ~{whitelist} > WhiteList.txt


    fastq_slideseq  \
      --bam-size 30.0 \
      --white-list WhiteList.txt \
      --read-structure "~{read_structure}" \
      --sample-id "~{sample_id}" \
      --output-format FASTQ \
      $FASTQS
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mb} MiB"
    disks: "local-disk ${disk} HDD"
    preemptible: preemptible
  }
  
  output {
    Array[File] fastq_R1_output_array = glob("fastq_R1_*")
    Array[File] fastq_R2_output_array = glob("fastq_R2_*")
  }
}

