version 1.0

task create_metadata_tsv {
    input {
        Array[File] sample_loom
        Array[String] library
        Array[String] species
        Array[String] stage
        Array[String] organ
        String output_tsv
        String docker = "quay.io/humancellatlas/secondary-analysis-loom-output:0.0.4-metadata-processing"
    }
    command <<<
    ~{sep='\n' sample_loom} >> sample_loom_file
    ~{sep='\n' library} >> library_file
    ~{sep='\n' species} >> species_file
    ~{sep='\n' stage} >> stage_file
    ~{sep='\n' organ} >> organ_file
    paste -d '\t' sample_loom_file library_file species_file stage_file organ_file > ~{output_tsv}
    >>>

    runtime {
        docker: docker
        cpu: 4  # note that only 1 thread is supported by pseudobam
        memory: "3 GiB"
        disks: "local-disk 100 HDD"
      }

    output {
    File metadata_tsv = "~{output_tsv}.tsv"
    }
}

workflow PostProcessingLoom {
  input {
    #runtime values
    Array[File] sample_loom
    Array[String] library
    Array[String] species
    Array[String] stage
    Array[String] organ
    String output_tsv = "metadata.tsv"

  }

  meta {
    description: "This task will add metadata to each loom file processed by our pipelines. This task will also generate json files for each loom file."
  }

  call create_metadata_tsv {
      input:
          sample_loom = sample_loom,
          library = library,
          species = species,
          stage = stage,
          organ = organ,
          output_tsv = output_tsv
  }

  output {
    File output_tsv = create_metadata_tsv.metadata_tsv
    ##File combined_loom = "~{output_loom}.loom"
    ##Array[File] modified_loom = "~{}.loom"
    ##Array[File] output_json = "~{}.json"
  }
}
