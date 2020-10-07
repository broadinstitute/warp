version 1.0

task FastqProcessing {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    File whitelist
    String chemistry
    String sample_id

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 3850
    Int cpu = 16   
    #TODO decided cpu
    # estimate that bam is approximately equal in size to fastq, add 20% buffer
    Int disk = ceil(size(r1_fastq, "GiB")*3 + size(r2_fastq, "GiB")*3) + 500

    Int preemptible = 3
  }

  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500

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

    # I1 file are optional,  and sometimes they are left out
    if [ -n '${sep=',' i1_fastq}' ]; then
      FLAG="--I1 ${sep=' --I1 ' i1_fastq}"
    else
      FLAG=''
    fi

    # use the right UMI length depending on the chemistry
    if [ "${chemistry}" == "tenX_v2" ]; then
        ## V2
        UMILENGTH=10
    elif [ "${chemistry}" == "tenX_v3" ]; then
        ## V3
        UMILENGTH=12
    else
        echo Error: unknown chemistry value: "$chemistry"
        exit 1;
    fi

    fastqprocess \
        --bam-size 3.0 \
        --barcode-length 16 \
        --umi-length $UMILENGTH \
        --sample-id "${sample_id}" \
        $FLAG \
        --R1 ${sep=' --R1 ' r1_fastq} \
        --R2 ${sep=' --R2 ' r2_fastq} \
        --white-list "${whitelist}"
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    Array[File] bam_output_array = glob("subfile_*")
  }
}
