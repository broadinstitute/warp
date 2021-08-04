version 1.0

# Get Cromwell metadata for the workflow that produced the given output
task get_metadata {
  input {
    String analysis_output_path
    String cromwell_url
    String include_keys = ""
    String include_subworkflows
    String pipeline_tools_version
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10

    String include_keys_arg = if (include_keys != "") then "--include_keys " + include_keys else ""
  }

  command <<<
    # Force the binary layer of the stdout and stderr streams (which is available as their buffer attribute)
    # to be unbuffered. This is the same as "-u", more info: https://docs.python.org/3/using/cmdline.html#cmdoption-u
    export PYTHONUNBUFFERED=TRUE

    get-analysis-workflow-metadata \
      --analysis_output_path ~{analysis_output_path} \
      --cromwell_url ~{cromwell_url} \
      --include_subworkflows ~{include_subworkflows} \
      ~{include_keys_arg}
  >>>
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-pipeline-tools:" + pipeline_tools_version
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"

  }
  output {
    File metadata = "metadata.json"
    String workflow_id = read_string("workflow_id.txt")
  }
}


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


task get_single_sample_ss2_inputs_from_metadata {
  input {
    File metadata_json
    String pipeline_tools_version
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command <<<
    # Force the binary layer of the stdout and stderr streams to be unbuffered.
    python -u <<CODE
    from pipeline_tools.pipelines.smartseq2 import smartseq2

    smartseq2.create_single_sample_ss2_inputs_tsv_from_analysis_metadata("~{metadata_json}")

    CODE
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-pipeline-tools:" + pipeline_tools_version
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    File inputs = "inputs.tsv"
  }

}

task get_inputs_from_metadata {
  input {
    File metadata_json
    String pipeline_tools_version
    String method
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  String import_statement = if method == 'MultiSampleSmartSeq2' then 'from pipeline_tools.pipelines.smartseq2 import smartseq2' else 'from pipeline_tools.pipelines.optimus import optimus'
  String run_statement = if method == 'MultiSampleSmartSeq2' then 'smartseq2.create_multi_sample_ss2_inputs_tsv_from_analysis_metadata' else 'optimus.create_optimus_inputs_tsv_from_analysis_metadata'
  command <<<
    # Force the binary layer of the stdout and stderr streams to be unbuffered.
    python -u <<CODE
    ~{import_statement}

    ~{run_statement}("~{metadata_json}")

    CODE
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-pipeline-tools:" + pipeline_tools_version
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    File inputs = "inputs.tsv"
    File? fastq1_input_files_tsv = "fastq1_input_files.tsv"
    File? fastq2_input_files_tsv = "fastq2_input_files.tsv"
    File? input_ids_tsv = "input_ids.tsv"

  }

}

# Create the analysis metadata
task create_submission {
  input {
    String workflow_id
    String pipeline_version
    String input_uuid
    File metadata_json
    String run_type
    String method
    String schema_url
    String analysis_process_schema_version
    String analysis_protocol_schema_version
    String analysis_file_version
    File inputs
    File outputs
    String pipeline_tools_version
    Boolean add_md5s
    String version
    String reference_uuid
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
    File? fastq1_input_files_tsv
    File? fastq2_input_files_tsv
    File? input_ids_tsv
  }

  command <<<
    # Force the binary layer of the stdout and stderr streams (which is available as their buffer attribute)
    # to be unbuffered. This is the same as "-u", more info: https://docs.python.org/3/using/cmdline.html#cmdoption-u
    export PYTHONUNBUFFERED=TRUE

    # First, create both analysis_process.json and analysis_protocol.json
    # Note that create-analysis-metadata can take a comma-separated list of bundles,
    # but current workflows only take a single input bundle
    create-analysis-metadata \
      --input_uuid "~{input_uuid}" \
      --analysis_id ~{workflow_id} \
      --metadata_json ~{metadata_json} \
      --run_type ~{run_type} \
      --method ~{method} \
      --schema_url ~{schema_url} \
      --analysis_process_schema_version ~{analysis_process_schema_version} \
      --analysis_protocol_schema_version ~{analysis_protocol_schema_version} \
      --pipeline_version ~{pipeline_version} \
      --analysis_file_version ~{analysis_file_version} \
      --inputs_file ~{inputs} \
      --outputs_file ~{outputs} \
      --add_md5s ~{add_md5s} \
      --version ~{version} \
      --references "~{reference_uuid}" \
      ~{"--fastq1_input_files_tsv=" + fastq1_input_files_tsv} \
      ~{"--fastq2_input_files_tsv=" + fastq2_input_files_tsv} \
      ~{"--input_ids_tsv=" + input_ids_tsv}
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-pipeline-tools:" + pipeline_tools_version
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    Array[File] analysis_process = glob("analysis_process/*.json")
    Array[File] analysis_protocol = glob("analysis_protocol/*.json")
    Array[File] analysis_files = glob("analysis_files/*.json")
    File outputs_file = "outputs.json"
    File inputs_file = "inputs.json"
  }
}


# Create the links json
task create_links {
  input{
    String links_schema_version
    String schema_url
    File analysis_process
    File analysis_protocol
    File outputs
    Array[String] input_uuids
    String pipeline_tools_version
    String version
    String project_id
    String project_stratum_string
    String input_id
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command <<<
    # Force the binary layer of the stdout and stderr streams (which is available as their buffer attribute)
    # to be unbuffered. This is the same as "-u", more info: https://docs.python.org/3/using/cmdline.html#cmdoption-u
    export PYTHONUNBUFFERED=TRUE

    # Build the links json
    create-links \
      --analysis_process_path ~{analysis_process} \
      --analysis_protocol_path ~{analysis_protocol} \
      --input_uuids ~{sep=' ' input_uuids} \
      --outputs_file_path ~{outputs} \
      --schema_url ~{schema_url} \
      --links_schema_version ~{links_schema_version} \
      --project_id ~{project_id} \
      --project_stratum_string "~{project_stratum_string}" \
      --input_id "~{input_id}" \
      --version ~{version}
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-pipeline-tools:" + pipeline_tools_version
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }
  output {
    Array[File] links_files = glob("links/*.json")
  }
}


task copy_file {
  input {
    String dest
    File src
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int additional_disk = 30
  }

  # 'gsutil cp' accepts either a file path or a "directory" path as its 2nd argument.
  # Since there's no real thing as a cloud directory, the only way it can tell which thing you mean is
  # if the destination ends with '/' or not.
  #
  # The `sub` call here will remove everything up to and including the final '/' in the input destination path.
  # If the result is an empty string, it means the destination ends with '/' and we're copying into a cloud
  # "directory". We could pass this directly into the `gsutil` call, but it's useful for downstream tasks
  # if we output the full destination path from this task, so we pre-compute a full file path to use.
  String dest_file = if sub(dest, "gs://.*/", "") == "" then dest + basename(src) else dest
  Int disk_size = ceil(size(src, "GiB")) + additional_disk
  command <<<
    mkdir gsutil_tmp
    export CLOUDSDK_CONFIG=gsutil_tmp
    gsutil -m cp ~{src} ~{dest_file}
  >>>

  output {
    Boolean done = true
    String copied_file = dest_file
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
  }
}


task get_cloud_file_creation_date {
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

task get_sha256_tsv {
  input {
    Array[String] file_paths
    Array[File] files
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int additional_disk = 30
  }
  Int disk_size = ceil(size(files[0], "GiB")) * length(files) + additional_disk

  command <<<
    files_array=(~{sep=' ' files})
    file_paths_array=(~{sep=' ' file_paths})
    for i in ${!files_array[*]}; do
      sha_output=$(sha256sum ${files_array[$i]})
      sha=$(echo $sha_output | awk '{print $1}')
      filename=$(echo $sha_output | awk '{print $2}')
      timestamp=$(gsutil ls -l ${file_paths_array[$i]} | egrep -o "([0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z)")
      echo $sha $filename $timestamp >> sha.tsv
    done
  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
  }

  output {
    File sha_tsv = "sha.tsv"
  }
}


task create_descriptor {
  input {
    String input_uuid
    File input_file
    String entity_type
    String input_file_path
    String schema_url
    String schema_version
    String pipeline_tools_version
    String creation_date
    String workspace_version
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int additional_disk = 30
  }

  Int disk_size = ceil(size(input_file, "GiB")) + additional_disk

  # TODO: Figure out how to get the actual file updated time within Terra
  command <<<
    export schema_url=~{schema_url}
    export schema_version=~{schema_version}
    export sha=$(sha256sum ~{input_file} | cut -f1 -d ' ')
    export crc=$(gsutil hash -h ~{input_file_path} | awk '/crc32c/ { print $3 }')
    export size=$(gsutil stat ~{input_file_path} | awk '/Content-Length/ { print $2 }')

    echo $sha > sha256.txt

    create-file-descriptor \
      --input_uuid "~{input_uuid}" \
      --file_path ~{input_file} \
      --entity_type ~{entity_type} \
      --size $size \
      --sha256 $sha \
      --crc32c $crc \
      --creation_time ~{creation_date} \
      --workspace_version ~{workspace_version} \
      --schema_url ~{schema_url} \
      --file_descriptor_schema_version ~{schema_version}
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-pipeline-tools:" + pipeline_tools_version
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
  }

  output {
    Array[File] descriptor_json = glob("*.json")
    String sha256 = read_string("sha256.txt")
  }
}


task create_reference_file {
  input {
    String file_path
    String input_uuid
    String schema_url
    String reference_file_schema_version
    String ncbi_taxon_id
    String version
    String reference_type
    String assembly_type
    String genus_species
    String reference_version
    String pipeline_tools_version
    Int cpu = 1
    Int machine_mem_mb = 2000
    Int disk = 10
  }

  command <<<
    create-reference-file \
      --file_path ~{file_path} \
      --input_uuid ~{input_uuid} \
      --schema_url ~{schema_url} \
      --reference_file_schema_version ~{reference_file_schema_version} \
      --ncbi_taxon_id ~{ncbi_taxon_id} \
      --workspace_version ~{version} \
      --reference_type '~{reference_type}' \
      --assembly_type '~{assembly_type}' \
      --genus_species '~{genus_species}' \
      --reference_version '~{reference_version}'
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-pipeline-tools:" + pipeline_tools_version
    cpu: cpu
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
  }

  output {
    Array[File] reference_file_json = glob("*.json")
    String reference_uuid = read_string("reference_uuid.txt")
  }
}


workflow submit {
  input{
    # Outputs from analysis
    Array[String] multisample_bams = []
    Array[String] multisample_bais = []
    Array[String] outputs

    # Constants shared across all files
    String run_type = "run"
    String schema_url = "https://schema.humancellatlas.org/"
    String cromwell_url = "https://api.firecloud.org/"

    # Schema versions for all generated metadata files
    String analysis_process_schema_version
    String analysis_protocol_schema_version
    String reference_file_schema_version
    String file_descriptor_schema_version
    String analysis_file_version
    String links_schema_version

    # Version of pipeline_tools to pull docker image
    String pipeline_tools_version

    Boolean add_md5s

    # Version of the pipeline, should be included in the pipeline file
    String method
    String pipeline_version
    String pipeline_version_string = method + "_v" + pipeline_version

    # Timestamp shared across analysis files in workspace
    String version = "2021-05-24T12:00:00.000000Z"

    # Project UUID
    Array[String] projects
    String project_id = projects[0]

    # Populate manually just to see if this works
    String genome_ref_fasta

    # More manual population for reference files
    String ncbi_taxon_id
    String genus_species
    String reference_type
    String assembly_type
    String reference_version

    # For copying files around
    String staging_area = "gs://broad-dsp-monster-hca-prod-lantern/"

    Array[String] input_uuids
    Array[String] fastq_1_uuids
    Array[String] fastq_2_uuids
    Array[String] fastq_i1_uuids = []

    Array[String] fastq_uuids = flatten([fastq_1_uuids, fastq_2_uuids, fastq_i1_uuids])
    Array[Array[String]] output_rows = if method == 'MultiSampleSmartSeq2' then transpose([fastq_1_uuids, fastq_2_uuids, multisample_bams, multisample_bais, input_uuids]) else []

    # Important metadata for UUID generation
    Array[String] library
    Array[String] species
    Array[String] organ

    # Use this to generate unique project-level UUIDs
    String project_stratum_string = "project=" + project_id + ";library=" + library[0] + ";species=" + species[0] + ";organ=" + organ[0]

    # Entity id source depends on pipeline
    String entity_project_identifier = if method == 'MultiSampleSmartSeq2' then project_stratum_string else input_uuids[0]

    # Build staging bucket
    String staging_bucket = staging_area + project_id + "/staging/"
  }

  call CheckMetadata {
    input:
      library = library,
      species = species,
      organ = organ
  }

  scatter (output_file in outputs) {
    call get_cloud_file_creation_date as get_creation_date_analysis_file_outputs {
      input:
        file_path = output_file
    }

    call create_descriptor as create_descriptor_analysis_file {
      input:
        input_uuid = entity_project_identifier,
        input_file = output_file,
        input_file_path = output_file,
        schema_url = schema_url,
        schema_version = file_descriptor_schema_version,
        pipeline_tools_version = pipeline_tools_version,
        creation_date = get_creation_date_analysis_file_outputs.creation_date,
        entity_type = "analysis_file",
        workspace_version = version
    }

    call copy_file as copy_file_analysis_output {
      input:
        dest = staging_bucket + "data/",
        src = output_file,
    }

    call copy_file as copy_file_analysis_file_descriptor {
      input:
        dest = staging_bucket + "descriptors/analysis_file/",
        src = create_descriptor_analysis_file.descriptor_json[0],
    }

  }

  scatter (output_row in output_rows) {
    call get_metadata as get_metadata_scatter {
      input:
        analysis_output_path = output_row[2],
        cromwell_url = cromwell_url,
        pipeline_tools_version = pipeline_tools_version,
        include_subworkflows = "True",
        include_keys = "input output calls start end runtimeAttributes stderr stdout"
    }

    call get_single_sample_ss2_inputs_from_metadata as get_inputs_from_ss2_bams {
      input:
        metadata_json = get_metadata_scatter.metadata,
        pipeline_tools_version = pipeline_tools_version
    }

    call create_submission as create_submission_ss2 {
      input:
        input_uuid = output_row[4],
        inputs = get_inputs_from_ss2_bams.inputs,
        run_type = run_type,
        schema_url = schema_url,
        analysis_process_schema_version = analysis_process_schema_version,
        analysis_protocol_schema_version = analysis_protocol_schema_version,
        analysis_file_version = analysis_file_version,
        method = method,
        outputs = get_sha256_tsv_ss2.sha_tsv,
        metadata_json = get_metadata_scatter.metadata,
        workflow_id = get_metadata_scatter.workflow_id,
        pipeline_version = pipeline_version,
        pipeline_tools_version = pipeline_tools_version,
        add_md5s = add_md5s,
        version = version,
        reference_uuid = create_reference_file_reference_genome.reference_uuid
    }

    call create_links as create_links_ss2 {
      input:
        schema_url = schema_url,
        analysis_process = create_submission_ss2.analysis_process[0],
        analysis_protocol = create_submission_ss2.analysis_protocol[0],
        input_uuids = [output_row[0], output_row[1]],
        outputs = create_submission_ss2.outputs_file,
        pipeline_tools_version = pipeline_tools_version,
        version = version,
        project_id = project_id,
        project_stratum_string = entity_project_identifier,
        input_id = output_row[4],
        links_schema_version = links_schema_version
    }

    call get_sha256_tsv as get_sha256_tsv_ss2 {
      input:
        file_paths = [output_row[2], output_row[3]],
        files = [output_row[2], output_row[3]],
    }

    scatter (analysis_file in create_submission_ss2.analysis_files) {
      call copy_file as copy_file_analysis_file_ss2 {
        input:
          dest = staging_bucket + "metadata/analysis_file/",
          src = analysis_file,
      }
    }

    call get_cloud_file_creation_date as get_creation_date_analysis_file_outputs_ss2_bam {
      input:
        file_path = output_row[2]
    }

    call get_cloud_file_creation_date as get_creation_date_analysis_file_outputs_ss2_bai {
      input:
        file_path = output_row[3]
    }

    call create_descriptor as create_descriptor_analysis_file_ss2_bam {
      input:
        input_uuid = output_row[4],
        input_file = output_row[2],
        input_file_path = output_row[2],
        schema_url = schema_url,
        schema_version = file_descriptor_schema_version,
        pipeline_tools_version = pipeline_tools_version,
        creation_date = get_creation_date_analysis_file_outputs_ss2_bam.creation_date,
        entity_type = "analysis_file",
        workspace_version = version
    }

    call create_descriptor as create_descriptor_analysis_file_ss2_bai {
      input:
        input_uuid = output_row[4],
        input_file = output_row[3],
        input_file_path = output_row[3],
        schema_url = schema_url,
        schema_version = file_descriptor_schema_version,
        pipeline_tools_version = pipeline_tools_version,
        creation_date = get_creation_date_analysis_file_outputs_ss2_bai.creation_date,
        entity_type = "analysis_file",
        workspace_version = version
    }

    call copy_file as copy_file_analysis_output_ss2_bam {
      input:
        dest = staging_bucket + "data/",
        src = output_row[2],
    }

    call copy_file as copy_file_analysis_output_ss2_bai {
      input:
        dest = staging_bucket + "data/",
        src = output_row[3],
    }

    call copy_file as copy_file_analysis_file_descriptor_ss2_bam {
      input:
        dest = staging_bucket + "descriptors/analysis_file/",
        src = create_descriptor_analysis_file_ss2_bam.descriptor_json[0],
    }

    call copy_file as copy_file_analysis_file_descriptor_ss2_bai {
      input:
        dest = staging_bucket + "descriptors/analysis_file/",
        src = create_descriptor_analysis_file_ss2_bai.descriptor_json[0],
    }

    call copy_file as copy_file_analysis_process_ss2 {
      input:
        dest = staging_bucket + "metadata/analysis_process/",
        src = create_submission_ss2.analysis_process[0],
    }

    call copy_file as copy_file_analysis_protocol_ss2 {
      input:
        dest = staging_bucket + "metadata/analysis_protocol/",
        src = create_submission_ss2.analysis_protocol[0],
    }

    call copy_file as copy_file_links_ss2 {
      input:
        dest = staging_bucket + "links/",
        src = create_links_ss2.links_files[0],
    }
  }

  call get_cloud_file_creation_date as get_creation_date_reference_file {
    input:
      file_path = genome_ref_fasta
  }

  call get_sha256_tsv {
    input:
      file_paths = outputs,
      files = outputs,
  }

  call get_metadata as get_metadata {
    input:
      analysis_output_path = outputs[0],
      cromwell_url = cromwell_url,
      pipeline_tools_version = pipeline_tools_version,
      include_subworkflows = "True",
      include_keys = "input output calls start end runtimeAttributes stderr stdout"

  }

  call get_metadata as get_metadata_inputs_only {
    input:
      analysis_output_path = outputs[0],
      cromwell_url = cromwell_url,
      pipeline_tools_version = pipeline_tools_version,
      include_subworkflows = "False",
      include_keys = "input"
  }


  call get_inputs_from_metadata {
    input:
      metadata_json = get_metadata.metadata,
      pipeline_tools_version = pipeline_tools_version,
      method = method
  }


  call create_descriptor as create_descriptor_reference_genome {
    input:
    # Reference files do not have an associated project, so we hash the filepath instead
      input_uuid = genome_ref_fasta,
      input_file = genome_ref_fasta,
      input_file_path = genome_ref_fasta,
      schema_url = schema_url,
      schema_version = file_descriptor_schema_version,
      pipeline_tools_version = pipeline_tools_version,
      creation_date = get_creation_date_reference_file.creation_date,
      entity_type = "reference_file",
      workspace_version = version
  }

  call create_reference_file as create_reference_file_reference_genome {
    input:
      file_path = genome_ref_fasta,
      input_uuid = genome_ref_fasta,
      schema_url = schema_url,
      reference_file_schema_version = reference_file_schema_version,
      ncbi_taxon_id = ncbi_taxon_id,
      version = version,
      reference_type = reference_type,
      assembly_type = assembly_type,
      genus_species = genus_species,
      reference_version = reference_version,
      pipeline_tools_version = pipeline_tools_version,
  }

  call create_submission {
    input:
      input_uuid = entity_project_identifier,
      inputs = get_inputs_from_metadata.inputs,
      run_type = run_type,
      schema_url = schema_url,
      analysis_process_schema_version = analysis_process_schema_version,
      analysis_protocol_schema_version = analysis_protocol_schema_version,
      analysis_file_version = analysis_file_version,
      method = method,
      outputs = get_sha256_tsv.sha_tsv,
      metadata_json = get_metadata.metadata,
      workflow_id = get_metadata.workflow_id,
      pipeline_version = pipeline_version,
      pipeline_tools_version = pipeline_tools_version,
      add_md5s = add_md5s,
      version = version,
      reference_uuid = create_reference_file_reference_genome.reference_uuid,
  }

  call create_links {
    input:
      schema_url = schema_url,
      analysis_process = create_submission.analysis_process[0],
      analysis_protocol = create_submission.analysis_protocol[0],
      outputs = create_submission.outputs_file,
      pipeline_tools_version = pipeline_tools_version,
      version = version,
      project_id = project_id,
      project_stratum_string = entity_project_identifier,
      input_id = "",
      links_schema_version = links_schema_version,
      input_uuids = fastq_uuids
  }

  call copy_file as copy_file_analysis_process {
    input:
      dest = staging_bucket + "metadata/analysis_process/",
      src = create_submission.analysis_process[0],
  }

  call copy_file as copy_file_analysis_protocol {
    input:
      dest = staging_bucket + "metadata/analysis_protocol/",
      src = create_submission.analysis_protocol[0],
  }

  scatter (analysis_file in create_submission.analysis_files) {
    call copy_file as copy_file_analysis_file {
      input:
        dest = staging_bucket + "metadata/analysis_file/",
        src = analysis_file,
    }
  }

  call copy_file as copy_file_reference_genome {
    input:
      dest = staging_bucket + "data/",
      src = genome_ref_fasta,
  }

  call copy_file as copy_file_reference_genome_reference_file {
    input:
      dest = staging_bucket + "metadata/reference_file/",
      src = create_reference_file_reference_genome.reference_file_json[0],
  }

  call copy_file as copy_file_reference_genome_descriptor {
    input:
      dest = staging_bucket + "descriptors/reference_file/",
      src = create_descriptor_reference_genome.descriptor_json[0],
  }

  call copy_file as copy_file_links {
    input:
      dest = staging_bucket + "links/",
      src = create_links.links_files[0],
  }


  output {
    Array[File] analysis_file = copy_file_analysis_file.copied_file
    File analysis_process = copy_file_analysis_process.copied_file
    File analysis_protocol = copy_file_analysis_protocol.copied_file
    Array[File] analysis_output = copy_file_analysis_output.copied_file
    File reference_genome = copy_file_reference_genome.copied_file
    File reference_genome_reference_file = copy_file_reference_genome_reference_file.copied_file
    File reference_genome_descriptor = copy_file_reference_genome_descriptor.copied_file
    Array[File] analysis_file_descriptor = copy_file_analysis_file_descriptor.copied_file
    File links = copy_file_links.copied_file

    # SS2 scatter outputs
    Array[File] descriptor_analysis_bam_ss2 = copy_file_analysis_file_descriptor_ss2_bam.copied_file
    Array[File] descriptor_analysis_bai_ss2 = copy_file_analysis_file_descriptor_ss2_bai.copied_file
    File descriptor_analysis_reference_ss2 = copy_file_reference_genome_descriptor.copied_file
    Array[File] analysis_file_ss2 = flatten(copy_file_analysis_file_ss2.copied_file)
    Array[File] analysis_process_ss2 = copy_file_analysis_process_ss2.copied_file
    Array[File] analysis_protocol_ss2 = copy_file_analysis_protocol_ss2.copied_file
    Array[File] links_ss2 = copy_file_links_ss2.copied_file
  }
}
