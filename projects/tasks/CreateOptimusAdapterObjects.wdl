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
    String creation_time
    String pipeline_type = "Optimus"
    String cromwell_url
  }

  call Tasks.GetCromwellMetadata {
    input:
      output_path = loom,
      cromwell_url = cromwell_url,
      include_subworkflows = false # TODO: do we need subworkflows???
  }

  # TODO: this will almost certainly fail for the projcet level run
  call Tasks.ParseCromwellMetadata {
    input:
      cromwell_metadata = GetCromwellMetadata.metadata,
      pipeline_type = pipeline_type
  }

  call Tasks.GetAnalysisFileMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      workspace_version = version_timestamp,
      input_file = GetCromwellMetadata.metadata
  }

  call Tasks.GetAnalysisProcessMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      workspace_version = version_timestamp,
      references = [],
      input_file = GetCromwellMetadata.metadata
  }

  call Tasks.GetAnalysisProtocolMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      workspace_version = version_timestamp,
      pipeline_version = ParseCromwellMetadata.pipeline_version
  }

  call Tasks.GetCloudFileCreationDate  as GetLoomFileCreationDate{
    input:
      file_path = loom
  }

  call Tasks.GetFileDescriptor as GetLoomFileDescriptor{
    input:
      pipeline_type = pipeline_type,
      file_path = loom,
      input_uuid = input_id,
      creation_time = GetLoomFileCreationDate.creation_date,
      workspace_version = version_timestamp
  }

  if (defined(bam)){
    call Tasks.GetCloudFileCreationDate  as GetBamFileCreationDate{
      input:
        file_path = select_first([bam])
    }

    call Tasks.GetFileDescriptor as GetBamFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = lselect_first([bam]),
        input_uuid = input_id,
        creation_time = GetBamFileCreationDate.creation_date,
        workspace_version = version_timestamp
    }
  }

  call Tasks.GetLinksFileMetadata {
    input:
      project_id = project_id,
      process_input_ids = process_input_ids, # for intermediate level use fastq_uuids from Terra, for project level use output_ids from intermediate files
      output_file_path = GetAnalysisFileMetadata.outputs_json,
      workspace_version = version_timestamp,
      analysis_process_path = GetAnalysisProcessMetadata.outputs_json,
      analysis_protocol_path = GetAnalysisProtocolMetadata.outputs_json,
      file_name_string = input_id
  }







#  call Tasks.GetDescriptorsAnalysisFileMetadata as GetDescriptorsAnalysisFileMetadataIntermediateLevelBam {
#    input:
#      pipeline_type = "Optimus",
#      file_path = bam,
#      input_uuid = input_id,
#      creation_time = GetCloudFileCreationDate.creation_date, #look at task get_cloud_file_creation_date in old wdl
#      workspace_version = version_timestamp
#  }
#
#  call Tasks.GetAnalysisFileMetadata as GetAnalysisFileMetadataProjectLevel {
#    input:
#      input_uuid = input_id,
#      pipeline_type = "Optimus",
#      workspace_version = version_timestamp,
#      input_file = MergeLooms.project_loom, #how do we get this
#      project_level = true
#  }
#
#  call Tasks.GetAnalysisProcessMetadata as GetAnalysisProcessMetadataProjectLevel {
#    input:
#      input_uuid = input_id,
#      pipeline_type = "Optimus",
#      workspace_version = version_timestamp,
#      references = [],
#      input_file = GetMetadata.metadata,
#      project_level=true,
#      loom_timestamp = #add this variable in TIMESTAMP=$(get_timestamp $LOOM_PATH)
#  }
#
#  call Tasks.GetAnalysisProtocolMetadata as GetAnalysisProtocolMetadataProjectLevel {
#    input:
#      input_uuid = input_id,
#      pipeline_type = "Optimus",
#      workspace_version = version_timestamp,
#      pipeline_version = GetMetadata.pipeline_version,
#      project_level = true
#  }
#
#  call Tasks.GetLinksFileMetadata as GetLinksFileMetadataProjectLevel {
#    input:
#      project_id = project_id,
#      process_input_ids = fastq_uuids, #come back to this, we want input_metadata.json for project level https://console.cloud.google.com/storage/browser/_details/fc-c307d7b3-8386-40a1-b32c-73b9e16e0103/b22deff9-924d-4aaa-a813-7a4d9d880915/TestHcaAdapter/215b754a-0657-45b8-a380-62db662b79a8/call-target_OptimusPostProcessing/OptimusPostProcessing/59d61014-c81f-4b05-b78e-395c62054a85/call-CreateAdapterJson/cacheCopy/script?authuser=0
#      output_file_path = GetAnalysisFileMetadataProjectLevel.outputs_json,
#      workspace_version = version_timestamp,
#      analysis_process_path = GetAnalysisProcessMetadataProjectLevel.outputs_json,
#      analysis_protocol_path = GetAnalysisProtocolMetadataProjectLevel.outputs_json,
#      file_name_string = project_stratum_string,
#      project_level=true
#  }

  output {
    File metadata_json = GetCromwellMetadata.metadata
    String reference_fasta = ParseCromwellMetadata.ref_fasta
    Array[File] analysis_file_outputs = GetAnalysisFileMetadata.analysis_file_outputs
    Array[File] analysis_process_outputs = GetAnalysisProcessMetadata.analysis_process_outputs
    Array[File] analysis_protocol_outputs = GetAnalysisProtocolMetadata.analysis_protocol_outputs
    Array[File] links_outputs = GetLinksFileMetadata.links_outputs
    Array[File] loom_file_descriptor_outputs = GetLoomFileDescriptor.file_descriptor_outputs
    Array[File]? bam_file_descriptor_outputs = GetBamFileDescriptor.file_descriptor_outputs
  }
}

