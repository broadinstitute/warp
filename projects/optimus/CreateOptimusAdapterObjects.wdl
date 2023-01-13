version 1.0

import "../tasks/AdapterTasks.wdl" as Tasks

workflow CreateOptimusAdapterObjects {
  meta {
    description: "Creates json objects for indexing HCA 10x analysis data"
    allowNestedInputs: true
  }

  input {
    File? bam
    File loom
    Array[String] process_input_ids # Array of space seperated strings...fastq for intermediate, intermediate looms for project level
    String input_id
    String project_id
    String version_timestamp
    String cromwell_url
    Boolean is_project_level
    String? pipeline_version # parsed from metadata for intermediate, passed in for project level
    String? reference_file_fasta # parsed from metadata for intermediate, passed in for project level
  }

  String pipeline_type = "Optimus"

  call Tasks.GetCromwellMetadata {
    input:
      output_path = loom,
      cromwell_url = cromwell_url,
      include_subworkflows = false
  }

  if (!is_project_level) {
    call Tasks.ParseCromwellMetadata {
      input:
        cromwell_metadata = GetCromwellMetadata.metadata,
        pipeline_type = pipeline_type
    }
  }

  String reference = select_first([ParseCromwellMetadata.ref_fasta, reference_file_fasta])
  String pipe_version = select_first([ParseCromwellMetadata.pipeline_version, pipeline_version])

  call Tasks.GetAnalysisFileMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      input_file = GetCromwellMetadata.metadata,
      project_level = is_project_level,
      project_loom = loom
  }

  call Tasks.GetAnalysisProcessMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      references = reference,
      project_level = is_project_level,
      input_file = GetCromwellMetadata.metadata
  }

  call Tasks.GetAnalysisProtocolMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      project_level = is_project_level,
      pipeline_version = pipe_version
  }

  call Tasks.GetCloudFileCreationDate as GetLoomFileCreationDate {
    input:
      file_path = loom
  }

  call Tasks.GetFileDescriptor as GetLoomFileDescriptor {
    input:
      pipeline_type = pipeline_type,
      file_path = loom,
      file_path_string = loom,
      input_uuid = input_id,
      creation_time = GetLoomFileCreationDate.creation_date,
      version_timestamp = version_timestamp
  }

  if (defined(bam)){
    call Tasks.GetCloudFileCreationDate  as GetBamFileCreationDate {
      input:
        file_path = select_first([bam])
    }

    call Tasks.GetFileDescriptor as GetBamFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = select_first([bam]),
        file_path_string = select_first([bam]),
        input_uuid = input_id,
        creation_time = GetBamFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }
  }

  call Tasks.GetLinksFileMetadata {
    input:
      project_id = project_id,
      process_input_ids = process_input_ids, # for intermediate level use fastq_uuids from Terra, for project level use output_ids from intermediate files
      output_file_path = GetAnalysisFileMetadata.outputs_json,
      version_timestamp = version_timestamp,
      analysis_process_path = GetAnalysisProcessMetadata.analysis_process_outputs,
      analysis_protocol_path = GetAnalysisProtocolMetadata.analysis_protocol_outputs,
      project_level = is_project_level,
      file_name_string = input_id,
      pipeline_type = pipeline_type
  }

  output {
    File metadata_json = GetCromwellMetadata.metadata
    String reference_fasta = reference
    String pipeline_version_string = pipe_version
    Array[File] analysis_file_outputs = GetAnalysisFileMetadata.analysis_file_outputs
    Array[File] analysis_process_outputs = GetAnalysisProcessMetadata.analysis_process_outputs
    Array[File] analysis_protocol_outputs = GetAnalysisProtocolMetadata.analysis_protocol_outputs
    Array[File] links_outputs = GetLinksFileMetadata.links_outputs
    Array[File] loom_file_descriptor_outputs = GetLoomFileDescriptor.file_descriptor_outputs
    Array[File]? bam_file_descriptor_outputs = GetBamFileDescriptor.file_descriptor_outputs
  }
}

