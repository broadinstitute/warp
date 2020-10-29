version 1.0

task create_metadata_tsv {
    input {
        Array[String] sample_loom
        Array[String] library
        Array[String] species
        Array[String] stage
        Array[String] organ
        String output_tsv
        String docker = "quay.io/humancellatlas/secondary-analysis-loom-output:0.0.4-metadata-processing"
    }
    command <<<
        echo "output_loom_file\tlibrary\tspecies\tstage\torgan" > ~{output_tsv}
        paste -d '\t' ~{write_lines(sample_loom)} ~{write_lines(library)} ~{write_lines(species)}  ~{write_lines(stage)} \
          ~{write_lines(organ)} >> ~{output_tsv}
    >>>

    runtime {
        docker: docker
        cpu: 1  # note that only 1 thread is supported by pseudobam
        memory: "3 GiB"
        disks: "local-disk 20 HDD"
      }

    output {
      File metadata_tsv = output_tsv
    }
}

workflow PostProcessingLoom {
  input {
    #runtime values
    Array[String] sample_loom
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
