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
    String library
    String species
    String input_id
    String organ
    String project_id
    String project_name
    String? project_stratum_string
    String version_timestamp
    String creation_time
    String genus_species #hoping to grab from metadata
    String reference_version #hoping to infer this, either GencodeV27 or GencodeM21
    String ncbi_taxon_id #hoping to infer this, either 9606 or 10090
    String assembly_type = "primary assembly" #hardcoded this for now
    String reference_type = "genome sequence" #hardcoded this for now
    String cromwell_url = "https://api.firecloud.org/"
  }

  call Tasks.GetCromwellMetadata as GetMetadata {
    input:
      output_path = loom,
      cromwell_url = cromwell_url,
      include_subworkflows = false # TODO: do we need subworkflows???
  }

  call Tasks.GetAnalysisFileMetadata as GetAnalysisFileMetadataIntermediateLevel {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      input_file = GetMetadata.metadata
  }
#need to parse metadata ahead of this step, make it available and we will pass it out and take it in here as references
  call Tasks.GetAnalysisProcessMetadata as GetAnalysisProcessMetadataIntermediateLevel {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      references = [],
      input_file = GetMetadata.metadata
  }

  call Tasks.GetAnalysisProtocolMetadata as GetAnalysisProtocolMetadataIntermediateLevel {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      pipeline_version = GetMetadata.pipeline_version
  }

  call Tasks.GetLinksFileMetadata as GetLinksFileMetadataIntermediateLevel {
    input:
      project_id = project_id,
      process_input_ids = fastq_uuids,
      output_file_path = GetAnalysisFileMetadata.outputs_json,
      workspace_version = version_timestamp,
      analysis_process_path = GetAnalysisProcessMetadata.outputs_json,
      analysis_protocol_path = GetAnalysisProtocolMetadata.outputs_json,
      file_name_string = input_id
  }

#  scatter (output_file in outputs) {
#    call Tasks.GetCloudFileCreationDate {
#      input:
#        file_path = output_files
#      }
#  }
  call Tasks.GetDescriptorsAnalysisFileMetadata as GetDescriptorsAnalysisFileMetadataIntermediateLevelLoom {
    input:
      pipeline_type = "Optimus",
      file_path = loom,
      input_uuid = input_id,
      creation_time = GetCloudFileCreationDate.creation_date, #look at task get_cloud_file_creation_date in old wdl
      workspace_version = version_timestamp
  }

  call Tasks.GetDescriptorsAnalysisFileMetadata as GetDescriptorsAnalysisFileMetadataIntermediateLevelBam {
    input:
      pipeline_type = "Optimus",
      file_path = bam,
      input_uuid = input_id,
      creation_time = GetCloudFileCreationDate.creation_date, #look at task get_cloud_file_creation_date in old wdl
      workspace_version = version_timestamp
  }

  call Tasks.GetAnalysisFileMetadata as GetAnalysisFileMetadataProjectLevel {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      input_file = MergeLooms.project_loom, #how do we get this
      project_level = true
  }

  call Tasks.GetAnalysisProcessMetadata as GetAnalysisProcessMetadataProjectLevel {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      references = [],
      input_file = GetMetadata.metadata,
      project_level=true,
      loom_timestamp = #add this variable in TIMESTAMP=$(get_timestamp $LOOM_PATH)
  }

  call Tasks.GetAnalysisProtocolMetadata as GetAnalysisProtocolMetadataProjectLevel {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      pipeline_version = GetMetadata.pipeline_version,
      project_level = true
  }

  call Tasks.GetLinksFileMetadata as GetLinksFileMetadataProjectLevel {
    input:
      project_id = project_id,
      process_input_ids = fastq_uuids, #come back to this, we want input_metadata.json for project level https://console.cloud.google.com/storage/browser/_details/fc-c307d7b3-8386-40a1-b32c-73b9e16e0103/b22deff9-924d-4aaa-a813-7a4d9d880915/TestHcaAdapter/215b754a-0657-45b8-a380-62db662b79a8/call-target_OptimusPostProcessing/OptimusPostProcessing/59d61014-c81f-4b05-b78e-395c62054a85/call-CreateAdapterJson/cacheCopy/script?authuser=0
      output_file_path = GetAnalysisFileMetadataProjectLevel.outputs_json,
      workspace_version = version_timestamp,
      analysis_process_path = GetAnalysisProcessMetadataProjectLevel.outputs_json,
      analysis_protocol_path = GetAnalysisProtocolMetadataProjectLevel.outputs_json,
      file_name_string = project_stratum_string,
      project_level=true
  }

  output {
    File metadata_json = GetMetadata.metadata
    Array[File] analysis_file_metadata = GetAnalysisFileMetadata.outputs
    Array[File] analysis_file_descriptor = GetDescriptorsAnalysisFileMetadata.ouptuts
    File analysis_process_metadata = GetAnalysisProcessMetadata.outputs
    File analysis_protocol_metdata = GetAnalysisProtocolMetadata.outputs
    File links = GetLinksFileMetadata.outputs
  }
}

