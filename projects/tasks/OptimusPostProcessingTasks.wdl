version 1.0

task CheckMetadata {
  input {
    Array[String] library
    Array[String] species
    Array[String] organ

    String docker = "python:3.7.2"
    Int memory = 3
    Int disk = 20
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

  if ';' in list(library_set)[0] or '=' in list(library_set)[0]:
      print('ERROR: Library metadata contains an illegal character (";" or "=")')
      errors += 1
  if ';' in list(species_set)[0] or '=' in list(species_set)[0]:
      print('ERROR: Species metadata contains an illegal character (";" or "=")')
      errors += 1
  if ';' in list(organ_set)[0] or '=' in list(organ_set)[0]:
      print('ERROR: Organ metadata contains an illegal character (";" or "=")')
      errors += 1

  if errors > 0:
      raise ValueError("Files must have matching metadata in order to combine.")
  CODE
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }
}


task MergeLooms {
  input {
    Array[File] library_looms
    String library
    String species
    String organ
    String project_id
    String project_name
    String output_basename

    String docker = "quay.io/humancellatlas/hca_post_processing:2.0"

    Int memory = ceil(size(library_looms, "G"))+ 10
    Int disk = ceil((size(library_looms, "G") * 4)) + 50
  }

  command {
    python3 /tools/optimus_HCA_loom_merge.py \
      --input-loom-files ~{sep=" " library_looms} \
      --library "~{library}" \
      --species "~{species}" \
      --organ "~{organ}" \
      --project-id "~{project_id}" \
      --project-name "~{project_name}" \
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


task GetInputMetadata {
  input {
    Array[File] analysis_file_jsons
    String output_basename

    String docker = "quay.io/humancellatlas/hca_post_processing:2.0"

    Int memory = ceil(size(analysis_file_jsons, "G")) + 3
    Int disk = ceil((size(analysis_file_jsons, "G") * 3)) + 20
  }
  command {
    python3 /tools/create_input_metadata_json.py \
      --input-json-files ~{sep=" " analysis_file_jsons} \
      --output ~{output_basename}.input_metadata.json
  }
  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    File input_metadata_json = "~{output_basename}.input_metadata.json"
  }
}


task CreateAdapterJson {
  input {
    File project_loom
    String project_id
    File input_metadata_json
    String project_stratum_string
    String staging_bucket
    String version_timestamp
    String pipeline_version

    Int memory = ceil((size(project_loom, "G") * 1.5)) + 5
    Int disk = ceil((size(project_loom, "G") * 3)) + 20
    String docker ="quay.io/humancellatlas/hca_post_processing:2.0"

  }

  command {
    source /tools/file_utils.sh

    LOOM_PATH=$(sed "s|/cromwell_root/|gs://|" <<< ~{project_loom})

    CRC=$(get_crc $LOOM_PATH)
    SHA=$(sha256sum ~{project_loom} | cut -f1 -d ' ')
    SIZE=$(get_size $LOOM_PATH)
    TIMESTAMP=$(get_timestamp $LOOM_PATH)

    mkdir outputs

    python3 /tools/HCA_create_adapter_json.py \
      --project-loom-file ~{project_loom} \
      --crc32c $CRC \
      --version-timestamp ~{version_timestamp} \
      --project-id ~{project_id} \
      --project-stratum-string "~{project_stratum_string}" \
      --sha256 $SHA \
      --size $SIZE \
      --staging-bucket ~{staging_bucket} \
      --input-metadata-json ~{input_metadata_json} \
      --loom-timestamp $TIMESTAMP \
      --pipeline-version ~{pipeline_version}
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    Array[File] json_adapter_files = glob("outputs/*")
  }
}
