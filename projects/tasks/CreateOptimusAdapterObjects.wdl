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

  call Tasks.GetAnalysisFileMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      metadata_json = GetMetadata.metadata
  }
#need to parse metadata ahead of this step, make it available and we will pass it out and take it in here as references
  call Tasks.GetAnalysisProcessMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = "Optimus",
      workspace_version = version_timestamp,
      references = ,
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
      pipeline_type = "Optimus",
      file_path = loom, #check to see if this runs for both loom or bam. if it just runs for bam then we would have to call it again
      input_uuid = input_id,
      creation_time = GetCloudFileCreationDate.creation_date, #look at task get_cloud_file_creation_date in old wdl
      workspace_version = version_timestamp
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

