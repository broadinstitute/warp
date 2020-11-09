version 1.0

task CheckMetadata {
  input {
    Array[String] library
    Array[String] species
    Array[String] organ
  }

  command {
  python <<CODE

  library_set = set([ "~{sep='", "' library}" ])
  species_set = set([ "~{sep='", "' species}" ])
  organ_set = set([ "~{sep='", "' organ}" ])

  errors=0

  if len(library_set) != 1:
      print("ERROR: Library metadata is not consistent within the project.")
      errors += 1
  if len(species_set) != 1:
      print("ERROR: Species metadata is not consistent within the project.")
      errors += 1
  if len(organ_set) != 1:
      print("ERROR: Organ metadata is not consistent within the project.")
      errors += 1

  if errors > 0:
      raise ValueError("Files must have matching metadata in order to combine.")
  CODE
  }

  runtime {
    docker: "python:3.7.2"
    cpu: 1
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }
}


task MergeLooms {
  input {
    Array[File] library_looms
    String library
    String species
    String organ
    String output_basename

    String docker = "quay.io/humancellatlas/secondary-analysis-loom-output:0.0.4-metadata-processing"
    Int memory = 3
    Int disk = 20
  }

  command {
    python3 optimus_HCA_loom_merge.py \
      --input-loom-files ~{sep=" " library_looms} \
      --library ~{library} \
      --species ~{species} \
      --organ ~{organ} \
      --output-loom-file ~{output_basename}.loom
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    File project_loom = "~{output_basename}.loom"
  }
}

task GetProtocolMetadata {
  input {
    Array[File] links_jsons
    String output_basename
  }
  command {
    python3 create_input_metadata_json.py \
      --input-files ~{sep=" " links_jsons} \
      --output ~{output_basename}.int_metadata.json
  }
  runtime {
    docker: "python:3.7.2"
    cpu: 1
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }
  output {
    File protocol_metadata_json = "~{output_basename}.protocol_metadata.json"
  }
}

task GetInputMetadata {
  input {
    Array[File] analysis_file_jsons
    String output_basename
  }
  command {
    python3 create_input_metadata_json.py \
      --input-files ~{sep=" " analysis_file_jsons} \
      --output ~{output_basename}.input_metadata.json
  }
  runtime {
    docker: "python:3.7.2"
    cpu: 1
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }
  output {
    File input_metadata_json = "~{output_basename}.input_metadata.json"
  }
}


task CreateAdapterJson {
  input {
    File project_loom
    File input_metadata_json
    String project_id

    Int memory = 3
    Int disk = 20
  }

  command {
    source file_utils.sh

    CRC=$(get_crc ~{project_loom})
    SHA=$(get_sha ~{project_loom})
    SIZE=$(get_size ~{project_loom})
    VERSION=$(get_timestamp ~{project_loom})

    python3 HCA_create_adapter_json.py \
      --input-loom-file ~{project_loom} \
      --inputs-json ~{input_metadata_json} \
      --project-id ~{project_id} \
      --crc32c $CRC \
      --size $SIZE \
      --sha256 $SHA \
      --version $VERSION \
  }

  runtime {
    docker: "python:3.7.2" # TODO make sure docker has gsutil
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    Array[File] json_adapter_files = glob("*$VERSION*.json")
  }
}