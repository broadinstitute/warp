version 1.0

task CheckMetadata {
  input {
    Array[String] library
    Array[String] species
    Array[String] stage
    Array[String] organ
  }

  #check that each peice of metadata is the same for all same for all indvidual looms
  command {
    LIBRARY=($(echo "~{library}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    SPECIES=($(echo "~{species}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    STAGE=($(echo "~{stage}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    ORGAN=($(echo "~{organ}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

    errors=0
    if [[ ${#LIBRARY[@]} != 1 ]]; then
      echo "ERROR: Library metadata is not consistent within the project."
      ((errors+=1))
    fi
    if [[ ${#SPECIES[@]} != 1 ]]; then
      echo "ERROR: Species metadata is not consistent within the project."
      ((errors+=1))
    fi
    if [[ ${#STAGE[@]} != 1 ]]; then
      echo "ERROR: Stage metadata is not consistent within the project."
      ((errors+=1))
    fi
    if [[ ${#ORGAN[@]} != 1 ]]; then
      echo "ERROR: Organ metadata is not consistent within the project."
      ((errors+=1))
    fi

    if [[ errors > 0 ]]; then
      exit 1
    fi
  }

  runtime { # TODO: change runtime setttings
      docker: docker
      cpu: 1  # note that only 1 thread is supported by pseudobam
      memory: "3 GiB"
      disks: "local-disk 20 HDD"
  }

  output {
    String library = "$LIBRARY"
    String species = "$SPECIES"
    String stage = "$STAGE"
    String organ = "$ORGAN"
  }

}


task CreateMetadataTsv {
  input {
    Array[String] original_looms
    Array[String] library
    Array[String] species
    Array[String] stage
    Array[String] organ
    String output_basename

    String docker = "quay.io/humancellatlas/secondary-analysis-loom-output:0.0.4-metadata-processing"
    Int memory = #TODO
    Int disk = #TODO
  }
  command <<<
      echo "output_loom_file  library species stage   organ" > ~{output_basename}.tsv
      paste -d '\t' ~{write_lines(sample_loom)} ~{write_lines(library)} ~{write_lines(species)}  ~{write_lines(stage)} \
        ~{write_lines(organ)} >> ~{output_basename}.tsv
  >>>

  runtime {
    docker: docker
    cpu: 1  # note that only 1 thread is supported by pseudobam
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }

  output {
    File metadata_tsv = "~{output_basename}.tsv"
  }
}


task IndividualLoomProcessing {
  input {
    File original_library_loom
    String library
    String species
    String stage
    String organ

    String docker = #TODO
    Int memory = #TODO
    Int disk = #TODO
  }
  String library_basename = basename(original_library_loom, ".loom")

  command {
    python3 optimus_individual_loom_post_processing.py \
      --input_loom ~{original_library_loom} \
      --library ~{library} \
      --species ~{species} \
      --stage ~{stage} \
      --organ ~{organ}
  }

  runtime { # TODO: change runtime setttings
    docker: docker
    cpu: 1  # note that only 1 thread is supported by pseudobam
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }

  output {
    File library_loom = "~{library_basename}.loom"
    File library_json = "~{library_basename}.json"
  }
}


task MergeLooms {
  input {
    Array[String] library_looms
    String output_basename

    String docker = "quay.io/humancellatlas/secondary-analysis-loom-output:0.0.4-metadata-processing"
    Int memory = #TODO
    Int disk = #TODO
  }

  command {
  }

  runtime { # TODO: change runtime setttings
    docker: docker
    cpu: 1  # note that only 1 thread is supported by pseudobam
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }

  output {
    File project_loom = "~{output_basename}.loom"
    File project_json = "~{output_basename}.json"
  }
}
