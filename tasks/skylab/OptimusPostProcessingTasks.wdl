version 1.0

task CheckMetadata {
  input {
    Array[String] library
    Array[String] species
    Array[String] organ
  }

  command {
  set -e pipefail
  python3 <<CODE

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
    String project_id

    String docker = "quay.io/humancellatlas/HCA_post_processing:0.0"
    Int memory = 3
    Int disk = 20
  }

  command {
    python3 tools/optimus_HCA_loom_merge.py \
      --input-loom-files ~{sep=" " library_looms} \
      --library ~{library} \
      --species ~{species} \
      --organ ~{organ} \
      --output-loom-file ~{project_id}.loom
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    File project_loom = "~{project_id}.loom"
  }
}


task GetInputMetadata {
  input {
    Array[File] analysis_file_jsons
    String project_id

    String docker = "quay.io/humancellatlas/HCA_post_processing:0.0"
  }
  command {
    python3 tools/create_input_metadata_json.py \
      --input-files ~{sep=" " analysis_file_jsons} \
      --output ~{project_id}.input_metadata.json
  }
  runtime {
    docker: docker
    cpu: 1
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }
  output {
    File input_metadata_json = "~{project_id}.input_metadata.json"
  }
}


task GetProtocolMetadata {
  input {
    Array[File] links_jsons
    String project_id

    String docker = "quay.io/humancellatlas/HCA_post_processing:0.0"
  }
  command {
    python3 tools/create_input_metadata_json.py \
      --input-files ~{sep=" " links_jsons} \
      --output ~{project_id}.protocol_metadata.json
  }
  runtime {
    docker: docker
    cpu: 1
    memory: "3 GiB"
    disks: "local-disk 20 HDD"
  }
  output {
    File protocol_metadata_json = "~{project_id}.protocol_metadata.json"
  }
}


task CreateAdapterJson {
  input {
    File project_loom
    String project_id
    File input_metadata_json
    File protocol_metadata_json
    String staging_bucket

    Int memory = 3
    Int disk = 20
    String docker ="quay.io/humancellatlas/HCA_post_processing:0.0"
  }

  command {
    source tools/file_utils.sh

    CRC=$(get_crc ~{project_loom})
    SHA=$(get_sha ~{project_loom})
    SIZE=$(get_size ~{project_loom})
    VERSION=$(get_timestamp ~{project_loom})

    python3 tools/HCA_create_adapter_json.py \
      --project-loom-file ~{project_loom} \
      --crc32c $CRC \
      --file_version $VERSION \
      --project-id ~{project_id} \
      --sha256 $SHA \
      --size $SIZE \
      --staging-bucket ~{staging_bucket} \
      --input-metadata-json ~{input_metadata_json} \
      --protocol-metadata-json ~{protocol_metadata_json}
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    Array[File] json_adapter_files = glob("*$VERSION*.json")
  }
}