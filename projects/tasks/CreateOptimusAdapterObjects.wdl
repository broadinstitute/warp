version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing

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

    String cromwell_url = "https://api.firecloud.org/"
  }

  call Tasks.GetCromwellMetadata as GetMetadata {
    input:
      output_path = loom,
      cromwell_url = cromwell_url,
      include_subworkflows = false # TODO: do we need subworkflows???
  }

  call Tasks.GetAnalysisFileMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      metadata_json = GetMetadata.metadata
  }

  call Tasks.GetAnalysisProcessMetadata {
    input:

  }

  call Tasks.GetAnalysisProtocolMetadata {
    input:

  }

  call Tsaks.GetLinksFileMetadata {
    input:
  }

  call Tsaks.GetDescriptorsAnalysisFileMetadata {
    input:
  }



  output {
    File metadata_json = GetMetadata.metadata
    Array[File] analysis_file_metadata = GetAnalysisFileMetadata.outputs
    Array[File] analysis_file_descriptor = GetDescriptorsAnalysisFileMetadata.ouptuts
    File analysis_process_metadata = GetAnalysisProcessMetadata.outputs
    File analysis_protocol_metdata = GetAnalysisProtocolMetadata.outputs
    File links = GetLinksFileMetadata.outtputs
  }
}

