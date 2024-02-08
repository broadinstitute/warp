version 1.0

workflow PIPseeker {

    input {
        String id #T1000-1
        Array[File] fastqs
        Array[File] snt_fastqs
        File snt_tags
        Array[File] refernce_files
    }

    # version of the pipeline
    String pipeline_version = "beta_0.0.1"

    call PIPseeker_full {
        input:
            id = id,
            fastqs = fastqs,
            snt_fastqs = snt_fastqs,
            snt_tags = snt_tags,
            refernce_files = refernce_files
    }

    output {
        Array[File] Results = PIPseeker_full.results
    }
}

task PIPseeker_full {
  input {
    String id
    Array[File] fastqs
    Array[File] snt_fastqs
    File snt_tags
    Array[File] refernce_files

    Int num_threads = 64
    String docker_image = "public.ecr.aws/w3e1n2j6/fluent-pipseeker:3.1.2"
    Int mem_size = ceil(size(fastqs, "GiB") + size(snt_fastqs, "GiB")) + (3 * num_threads) #1GB for every 1GB fastq input in gz format + 3GB of ram for every additional processing thread
    Int disk_size = ceil((size(fastqs, "GiB") + size(snt_fastqs, "GiB")) * 2.5) + (3 * num_threads) #2.5GB for every 1GB fastq input in gz format + uncompressed size of ref genome
    Int preemptible_tries = 3
  }

  command <<<
    set -euo pipefail

    # SETUP
    # ensure the the fastq files all begin with a common prefix
    mkdir FASTQS
    declare -a FASTQ_ARRAY=(~{sep=' ' fastqs})
    for f in "${FASTQ_ARRAY[@]}"; do mv $f FASTQS; done

    # ensure the the snt-fastq files all begin with a common prefix
    mkdir SNT_FASTQS
    declare -a SNT_FASTQ_ARRAY=(~{sep=' ' snt_fastqs})
    for f in "${SNT_FASTQ_ARRAY[@]}"; do mv $f SNT_FASTQS; done

    # ensure the the snt-fastq files all begin with a common prefix
    mkdir REFERNCE
    declare -a REF_ARRAY=(~{sep=' ' refernce_files})
    for f in "${REF_ARRAY[@]}"; do mv $f REFERNCE; done

    mkdir results

    # Run PIPseeker
    pipseeker-v3.1.2 full  \
      --fastq FASTQS/. \
      --id ~{id} \
      --verbosity 2 \
      --skip-preflight \
      --star-index-path REFERNCE \
      --threads ~{num_threads} \
      --monitor 10 \
      --chemistry v4  \
      --snt-fastq SNT_FASTQS/. \
      --snt-tags ~{snt_tags} \
      --output-path results
  >>>

  runtime {
    docker: docker_image
    disk: "1600 GiB" # "local-disk ${disk_size} HDD"
    cpu: 64 # cpu
    memory: "432 GiB"# "${mem_size} GiB"
    preemptible: preemptible_tries
    vm_size: "Standard_E64_v3"
  }

  output {
    Array[File] results = glob("results/*")
    }
}
