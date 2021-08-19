version 1.0

# Due to the nature generating adapters for intermediatt and project level files, we have several place
# where we process arrays of inputs which we expect to contain a single value duplicated many times.
# This simple function confirms that there is a sigle value in a given input array and that that value does not
# contain any disallowed characters
task CheckInput {
  input {
    Array[String] input_array
    String input_type
    String illegal_characters = ""

    String docker = "python:3.7.2"
    Int memory = 3
    Int disk = 10
  }

  command <<<
  set -e pipefail
  python3 <<CODE

  input_set = set([ "~{sep='", "' input_array}" ])

  errors=0

  if len(input_set) != 1:
      print("ERROR: Expected one value for {}, but found multiple: {}".format(input_type, input_set))
      errors += 1

  for c in list(~{illegal_characters}):
      for i in input_set:
          if c in i:
              print("ERROR: {} string, {}, contains an illegal character {}".format(input_type, i, c))
              errors += 1

  if errors > 0:
      raise ValueError("Failed input check due to mulitple values and/or illegal characters")

  with open('output.txt', 'w') as f:
        f.write(list(input_set)[0])
  CODE
  >>>
  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    String output_string = read_string("output.txt")
  }
}


task GetPipelineType {
  input {
    String library
    String docker = "python:3.7.2"
    Int memory = 3
    Int disk = 10
  }

  command <<<
  set -e pipefail
  python3 <<CODE
  with open("output.txt", w) as f:
      if ("10X" in ~{library}):
          f.write("Optimus")
      elif ("Smart-seq2" in ~{library}):
          f.write("SS2")
      else:
          raise ValueError("Unexpected library_preparation_protocol__library_construction_approach")
  CODE
  >>>
  runtime {
    docker: docker
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    String output_string = read_string("output.txt")
  }
}


# Get Cromwell metadata for a workflow
# Uses a workflow output to parse the cromwell id and fetch the metadata
task GetCromwellMetadata {
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


task MergeLooms {
  input {
    Array[File] output_looms
    String library
    String species
    String organ
    String project_id
    String project_name
    String output_basename

    String docker = "quay.io/humancellatlas/hca_post_processing:2.0"
    Int memory = ceil(size(output_looms, "G"))+ 10
    Int disk = ceil((size(output_looms, "G") * 4)) + 50
  }

  command {
    python3 /tools/optimus_HCA_loom_merge.py \
      --input-loom-files ~{sep=" " output_looms} \
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


task GetAnalysisFileMetadata {
  input {
    String input_uuid
    String pipeline_type
    String version_timestamp
    String input_file
    Boolean? project_level

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    create-analysis-file \
      --input_uuid = "~{input_uuid}" \
      --pipeline_type = "~{pipeline_type}" \
      --workspace_version = "~{version_timestamp}" \
      --input_file = "~{input_file}" \
      ~{"--project_level " + project_level}
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    Array[File] analysis_file_outputs = glob("*${version_timestamp}.json")
    File outputs_json = "outputs.json"
  }
}


task GetAnalysisProcessMetadata {
  input {
    String input_uuid
    String pipeline_type
    String version_timestamp
    String references
    String input_file
    Boolean? project_level
    String? loom_timestamp

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    create-analysis-process \
      --input_uuid = "~{input_uuid}" \
      --pipeline_type = "~{pipeline_type}" \
      --workspace_version = "~{version_timestamp}" \
      --references ="~{references}" \
      --input_file ="~{input_file}" \
      ~{"--project_level " + project_level} \
      ~{"--loom_timestamp " + loom_timestamp}

  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    Array[File] analysis_process_outputs = glob("*${version_timestamp}.json")
  }
}

task GetAnalysisProtocolMetadata {
   input {
     String input_uuid
     String pipeline_type
     String version_timestamp
     String pipeline_version
     Boolean? project_level

     String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
     Int cpu = 1
     Int machine_mem_mb = 2000
     Int disk = 10
   }

   command {
     create-analysis-protocol \
       --input_uuid = "~{input_uuid}" \
       --pipeline_type = "~{pipeline_type}" \
       --workspace_version = "~{version_timestamp}" \
       --pipeline_version = "~{pipeline_version}" \
       ~{"--project_level " + project_level}
   }
   runtime {
     docker: docker
     cpu: cpu
     memory: "${machine_mem_mb} MiB"
     disks: "local-disk ~{disk} HDD"
    }
    output {
      Array[File] analysis_protocol_outputs = glob("*${version_timestamp}.json")
    }
 }


task GetLinksFileMetadata {
  input {
    String project_id
    Array[String] process_input_ids
    String output_file_path
    String version_timestamp
    Array[String] analysis_process_path
    Array[String] analysis_protocol_path
    String file_name_string
    Boolean? project_level

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    create-links \
    --project_id = "~{project_id}" \
    --input_uuids = "~{sep=' ' process_input_ids}" \
    --output_file_path = "~{output_file_path}" \
    --workspace_version = "~{version_timestamp}" \
    --analysis_process_path = "~{sep=' ' analysis_process_path}" \
    --analysis_protocol_path = "~{sep=' ' analysis_protocol_path}" \
    --file_name_string = "~{file_name_string}" \
    ~{"--project_level " + project_level}
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    Array[File] links_outputs = glob("*${version_timestamp}*.json")
  }
}


task GetFileDescriptor {
  input {
    String input_uuid
    String pipeline_type
    String creation_time
    String version_timestamp
    File file_path
    String file_path_string #does this need to be set to file_path ?

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 30
   }

  command
  <<<
      export sha=$(sha256sum ~{file_path} | cut -f1 -d ' ')
      export crc=$(gsutil hash -h ~{file_path_string} | awk '/crc32c/ { print $3 }')
      export size=$(gsutil stat ~{file_path_string} | awk '/Content-Length/ { print $2 }')

    create-file-descriptor \
    --size = "$size" \
    --sha256 = "$sha256" \
    --crc32c = "$crc32c" \
    --pipeline_type = "~{pipeline_type}" \
    --file_path = "~{file_path}" \
    --input_uuid = "~{input_uuid}" \
    --creation_time = "~{creation_time}" \
    --workspace_version = "~{version_timestamp}"

  >>>
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    Array[File] file_descriptor_outputs = glob("*${version_timestamp}.json")
  }
}


task GetReferenceFileMetadata {
  input {
    String file_path
    String input_uuid
    String genus_species
    String assembly_type
    String pipeline_type
    String ncbi_taxon_id
    String reference_type
    String version_timestamp
    String reference_version

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
  create-reference-file \
  --genus_species = "~{genus_species}" \
  --file_path = "~{file_path}" \
  --workspace_version = "~{version_timestamp}" \
  --input_uuid = "~{input_uuid}" \
  --reference_version = "~{reference_type}" \
  --ncbi_taxon_id = "~{ncbi_taxon_id}" \
  --pipeline_type = "~{pipeline_type}" \
  --assembly_type = "~{assembly_type}" \
  --reference_type = "~{reference_type}"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    String reference_file_uuid = read_string("reference_uuid.txt")
    Array[File] reference_metadata_outputs = glob("*${version_timestamp}.json")
  }
}


task GetCloudFileCreationDate {
  input {
    String file_path
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command <<<
    gsutil ls -l ~{file_path} | egrep -o "([0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z)" > creation_date.txt
  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    String creation_date = read_string("creation_date.txt")
  }
}


task ParseCromwellMetadata {
  input {
    File cromwell_metadata
    String pipeline_type

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    parse-metadata \
    --cromwell-metadata-json ~{cromwell_metadata} \
    --pipeline-type ~{pipeline_type}
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    String ref_fasta = read_string("ref_fasta.txt")
    String pipeline_version = read_string("pipeline_version.txt")
  }
}


task GetReferenceDetails {
  input {
    File ref_fasta
    String species

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    get-reference-file-details \
    --reference-file ~{ref_fasta} \
    --species ~{species}
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    String ncbi_taxon_id = read_string("ncbi_taxon_id.txt")
    String assembly_type = read_string("assembly_type.txt")
    String reference_type = read_string("reference_type.txt")
  }
}

task GetProjectLevelInputIds {
  input {
    Array[File] intermediate_analysis_files

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    python3 get_process_input_ids.py \
    --input-json-files ~{sep=' ' intermediate_analysis_files}
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    String process_input_uuids = read_string("output.txt")
  }
}

task CopyToStagingBucket {
  input  {
    Array[File] analysis_file_metadata_objects
    Array[File] analysis_process_objects
    Array[File] analysis_protocol_objects
    Array[File] analysis_file_descriptor_objects
    Array[File] links_objects
    Array[File] data_objects
    String staging_bucket
    String? cache_invalidate

    String docker = "quay.io/humancellatlas/secondary-analysis-pipeline-tools:master"
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command {
    copy-adapter-outputs \
    --analysis_files_metadata_jsons ~{analysis_file_metadata_objects} \
    --analysis_process_jsons ~{analysis_process_objects} \
    --analysis_protocol_jsons ~{analysis_protocol_objects} \
    --analysis_files_descriptors_jsons ~{analysis_file_descriptor_objects} \
    --links_jsons ~{links_objects} \
    --data_files ~{data_objects} \
    --staging-bucket ~{staging_bucket}
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
  }
}

