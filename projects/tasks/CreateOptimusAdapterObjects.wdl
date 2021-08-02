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
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      references =,
      metadata_json = GetMetadata.metadata

  }

  call Tasks.GetAnalysisProtocolMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      pipeline_version = GetMetadata.pipeline_version
  }

  call Tasks.GetLinksFileMetadata {
    input:
      project_id = project_id,
      process_input_ids = fastq_uuids,
      output_file_path = GetAnalysisFileMetadata.outputs_json,
      workspace_version = version_timestamp,
      analysis_process_path = GetAnalysisProcessMetadata.outputs_json,
      analysis_protocol_path = GetAnalysisProtocolMetadata.outputs_json,
      file_name_string = input_id
  }

  call Tasks.GetDescriptorsAnalysisFileMetadata {
    input:
      size = ,
      sha256 = ,
      crc32c = ,
      pipeline_type = "Optimus",
      file_path = ,
      input_uuid = ,
      creation_time = ,
      workspace_version = version_timestamp
  }

  call Tasks.GetReferenceFileMetadata {
    input:
      genus_species = ,
      file_path = ,
      workspace_version = version_timestamp ,
      input_uuid = ,
      reference_version = ,
      ncbi_taxon_id = ,
      pipeline_type = "Optimus",
      assembly_type = ,
      reference_type =
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

