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

    for f in ~{sep=' ' r1_fastq}; do
      if [[ $f != *.fastq* ]]; then
        mv  "$f" "$f".fastq
        FQ="$f".fastq
      else
        FQ="$f"
      fi
      if [[ $FQ != *.gz ]]; then
        if (file $FQ | grep -q compressed); then
          mv  "$FQ" "$FQ".gz
          FQ="$FQ".gz
        fi
      fi
      R1_FASTQS+=" --R1 $FQ "
    done

    for f in ~{sep=' ' r2_fastq}; do
      if [[ $f != *.fastq* ]]; then
        mv  "$f" "$f".fastq
        FQ="$f".fastq
      else
        FQ="$f"
      fi
      if [[ $FQ != *.gz ]]; then
        if (file $FQ | grep -q compressed); then
          mv  "$FQ" "$FQ".gz
          FQ="$FQ".gz
        fi
      fi
      R2_FASTQS+=" --R2 $FQ "
    done

    # I1 file are optional,  and sometimes they are left out
    if [ -n '~{sep=',' i1_fastq}' ]; then
      for f in ~{sep=' ' i1_fastq}; do
        if [[ $f != *.fastq* ]]; then
          mv  "$f" "$f".fastq
          FQ="$f".fastq
        else
          FQ="$f"
        fi
        if [[ $FQ != *.gz ]]; then
          if (file $FQ | grep -q compressed); then
            mv  "$FQ" "$FQ".gz
            FQ="$FQ".gz
          fi
        fi
        I1_FASTQS+=" --I1 $FQ "
      done
    fi

    # use the right UMI length depending on the chemistry
    if [ "~{chemistry}" == "tenX_v2" ]; then
        ## V2
        UMILENGTH=10
    elif [ "~{chemistry}" == "tenX_v3" ]; then
        ## V3
        UMILENGTH=12
    else
        echo Error: unknown chemistry value: "~{chemistry}"
        exit 1;
    fi

    fastqprocess \
        --bam-size 1.0 \
        --barcode-length 16 \
        --umi-length $UMILENGTH \
        --sample-id "~{sample_id}" \
        "$I1_FASTQS" \
        "$R1_FASTQS" \
        "$R2_FASTQS" \
        --white-list "~{whitelist}"
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
