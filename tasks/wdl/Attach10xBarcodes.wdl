version 1.0

task Attach10xBarcodes {
  input {
    File r1_fastq
    File? i1_fastq
    File r2_unmapped_bam
    File whitelist
    String chemistry

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.11"
    Int machine_mem_mb = 48000
    Int cpu = 2
    # estimate that bam is approximately the size of all inputs plus 50%
    Int disk = ceil((size(r2_unmapped_bam, "Gi") + size(r1_fastq, "Gi") + if (defined(i1_fastq)) then size(i1_fastq, "Gi") else 0) * 3)
    # by default request non preemptible machine to make sure the slow attach barcodes step completes
    Int preemptible = 0
  }

  meta {
    description: "attaches barcodes found in r1 (forward) and i1 (index) fastq files to corresponding reads in the r2 (reverse) bam file"
  }

  parameter_meta {
    r1_fastq: "forward fastq file; contains umi, cell barcode"
    i1_fastq: "optional, index fastq file; contains sample barcode"
    r2_unmapped_bam: "reverse unmapped bam file; contains alignable genomic information"
    whitelist: "10x genomics cell barcode whitelist"
    chemistry: "chemistry employed, currently can be tenX_v2, tenX_v3, or v4 (GEM-X), the latter two implies NO feature barcodes"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    if [ "${chemistry}" == "tenX_v2" ]
    then
        ## V2
        Attach10xBarcodes \
            --r1 "${r1_fastq}" \
            ${"--i1 " + i1_fastq} \
            --u2 "${r2_unmapped_bam}" \
            --whitelist "${whitelist}" \
            --output-bamfile barcoded.bam
    elif [ "${chemistry}" == "tenX_v3" ] || [ "${chemistry}" == "tenX_v4" ]
    then
        ## V3
        AttachBarcodes \
            --r1 "${r1_fastq}" \
            ${"--i1 " + i1_fastq} \
            --u2 "${r2_unmapped_bam}" \
            --whitelist "${whitelist}" \
            --sample-barcode-start-position 0 \
            --sample-barcode-length 8 \
            --cell-barcode-start-position 0 \
            --cell-barcode-length 16 \
            --molecule-barcode-start-position 16 \
            --molecule-barcode-length 12 \
            --output-bamfile barcoded.bam
    else
        echo Error: unknown chemistry value: "$chemistry"
        exit 1;
    fi

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
    File bam_output = "barcoded.bam"
  }
}
