version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing
import "../../projects/tasks/AdapterTasks.wdl" as Tasks

workflow CreateOptimusAdapterObjects {
  meta {
    description: "Creates json objects for indexing HCA 10x analysis data"
    allowNestedInputs: true
  }

  input {
    File? bam
    File loom
    Array[String] process_input_ids # TODO: create array of strings for project level inputs (intermediate outputs)
    String library
    String species
    String input_id
    String organ
    String project_id
    String project_name
    String? project_stratum_string
    String version_timestamp
    String pipeline_type = "Optimus"
    String cromwell_url
    Boolean is_project_level
    String? pipeline_version # parsed from metadata for intermediate, passed in for project level
    String? reference_fasta # parsed from metadata for intermediate, passed in for project level
  }

  call Tasks.GetCromwellMetadata {
    input:
      output_path = loom,
      cromwell_url = cromwell_url,
      include_subworkflows = false # TODO: do we need subworkflows???
  }

  if (!is_project_level) {
    call Tasks.ParseCromwellMetadata {
      input:
        cromwell_metadata = GetCromwellMetadata.metadata,
        pipeline_type = pipeline_type
    }
  }
    String reference = select_first([ParseCromwellMetadata.ref_fasta,reference_fasta])

  call Tasks.GetAnalysisFileMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      input_file = GetCromwellMetadata.metadata
  }

  call Tasks.GetAnalysisProcessMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      references = select_first([ParseCromwellMetadata.ref_fasta,reference_fasta]),
      input_file = GetCromwellMetadata.metadata
  }

  call Tasks.GetAnalysisProtocolMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      pipeline_version = select_first([ParseCromwellMetadata.pipeline_version,pipeline_version])
  }

  call Tasks.GetCloudFileCreationDate  as GetLoomFileCreationDate {
    input:
      file_path = loom
  }

  call Tasks.GetFileDescriptor as GetLoomFileDescriptor {
    input:
      pipeline_type = pipeline_type,
      file_path = loom,
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
        input_uuid = input_id,
        creation_time = GetBamFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }
  }

  call Tasks.GetLinksFileMetadata {
    input:
      project_id = project_id,
      process_input_ids = process_input_ids, # for intermediate level use fastq_uuids from Terra, for project level use output_ids from intermediate files
      output_file_path = GetAnalysisFileMetadata.analysis_file_outputs,
      version_timestamp = version_timestamp,
      analysis_process_path = GetAnalysisProcessMetadata.analysis_process_outputs,
      analysis_protocol_path = GetAnalysisProtocolMetadata.analysis_protocol_outputs,
      file_name_string = input_id
  }


  output {
    File metadata_json = GetCromwellMetadata.metadata
    String reference_fasta = reference
    Array[File] analysis_file_outputs = GetAnalysisFileMetadata.analysis_file_outputs
    Array[File] analysis_process_outputs = GetAnalysisProcessMetadata.analysis_process_outputs
    Array[File] analysis_protocol_outputs = GetAnalysisProtocolMetadata.analysis_protocol_outputs
    Array[File] links_outputs = GetLinksFileMetadata.links_outputs
    Array[File] loom_file_descriptor_outputs = GetLoomFileDescriptor.file_descriptor_outputs
    Array[File]? bam_file_descriptor_outputs = GetBamFileDescriptor.file_descriptor_outputs
  }
}

