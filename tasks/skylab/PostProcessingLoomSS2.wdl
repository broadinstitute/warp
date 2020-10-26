version 1.0

workflow PostProcessingLoom {
  input {
    #runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-loom-output:0.0.3"
    File metadata_tsv
    String output_loom
  }

  meta {
    description: "This task will add metadata to each loom file processed by our pipelines. This task will also generate json files for each loom file."
  }

  command {
  #TODO: get first col of tsv file and pass as output_json name

    set -euo pipefail

    python3 /tools/post_processing_ss2.py \
       --input-file ~{metadata_tsv} \
       --output-loom-file ~{output_loom}
  }

  runtime {
    docker: docker
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "3 GiB"
    disks: "local-disk 100 HDD"
    preemptible: preemptible
  }

  output {
    File combined_loom = "~{output_loom}.loom"
    Array[File] modified_loom = "~{}.loom"
    Array[File] output_json = "~{}.json"
  }
}

