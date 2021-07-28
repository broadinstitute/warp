version 1.0

# Check that the library prep method, species, and organ are the same for all workflow being processed together
task CheckStratumMetadata {
  input {
    Array[String] library
    Array[String] species
    Array[String] organ
    Array[String] project_id
    Array[String] project_name
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
  project_id_set = set([ "~{sep='", "' project_id}" ])
  project_name_set = set([ "~{sep='", "' project_name}" ])

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
  if len(project_id_set) != 1:
        print("ERROR: Project_id metadata is not consistent within the project.")
        errors += 1
  if len(project_name_set) != 1:
      print("ERROR: Project_name metadata is not consistent within the project.")
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
  if ';' in list(project_id_set)[0] or '=' in list(project_id_set)[0]:
      print('ERROR: Project_id metadata contains an illegal character (";" or "=")')
      errors += 1
  if ';' in list(project_name_set)[0] or '=' in list(project_name_set)[0]:
      print('ERROR: Project_name metadata contains an illegal character (";" or "=")')
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

# Get Cromwell metadata for a workflow
# Uses a workflow output to parse the cromwell id and fetch the metadata
task GetMetadata {
  input {
    String output_path
    String cromwell_url
    Boolean include_subworkflows
    Array[String] include_keys

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:gw_cromwell_includekeys"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    # Force the binary layer of the stdout and stderr streams (which is available as their buffer attribute)
    # to be unbuffered. This is the same as "-u", more info: https://docs.python.org/3/using/cmdline.html#cmdoption-u
    export PYTHONUNBUFFERED=TRUE

    get-analysis-workflow-metadata \
      --analysis_output_path ~{output_path} \
      --cromwell_url ~{cromwell_url} \
      ~{true="--include_subworkflows True" false="--include_subworkflows False" include_subworkflows} \
      ~{"--include_keys " + include_keys}
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    File metadata = "metadata.json"
    String workflow_id = read_string("workflow_id.txt")
  }
}


