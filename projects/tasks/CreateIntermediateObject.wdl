version 1.0

# This workflow functions as a wrapper around the scatter call to CreateIntermediateOptimusAdapters
# As files are created in each iteration (each individual run of Optimus) we have to gather them all and export
# See the 'scatter-gather' pattern here https://github.com/openwdl/wdl/blob/main/versions/development/SPEC.md#advanced-wdl-features

import "../../projects/tasks/CreateOptimusAdapterObjects.wdl" as CreateOptimusObjects

workflow CreateIntermediateObject {
  input {
    Array[File] output_bams
    Array[File] output_looms
    Array[File]? output_bais
    Array[String] input_ids 
    Array[String] fastq_1_uuids
    Array[String] fastq_2_uuids
    Array[String]? fastq_i1_uuids

    String library 
    String organ 
    String species 
    String project_id 
    String project_name 

  
    String cromwell_url = "https://api.firecloud.org/"
    String version_timestamp = "2021-05-24T12:00:00.000000Z" 
  }

  if (false) {
    String none = "None"
  }

  scatter (idx in range(length(output_looms))) {
      String? fastq_i1_uuid = if defined(fastq_i1_uuids) then select_first([fastq_i1_uuids])[idx] else none
      call CreateOptimusObjects.CreateOptimusAdapterObjects as CreateIntermediateOptimusAdapters {
        input:
          bam = output_bams[idx],
          loom = output_looms[idx],
          input_id = input_ids[idx],
          process_input_ids = select_all([fastq_1_uuids[idx],fastq_2_uuids[idx], fastq_i1_uuid]),
          library = library,
          species = species,
          organ = organ,
          project_id = project_id,
          project_name = project_name,
          version_timestamp = version_timestamp,
          cromwell_url = cromwell_url,
          is_project_level = false
      }
  }

  output {
    Array[String] reference_fasta = CreateIntermediateOptimusAdapters.reference_fasta
    Array[String] pipeline_version_string = CreateIntermediateOptimusAdapters.pipeline_version_string
    Array[Array[File]] analysis_file_outputs = CreateIntermediateOptimusAdapters.analysis_file_outputs
    Array[Array[File]] analysis_process_outputs = CreateIntermediateOptimusAdapters.analysis_process_outputs
    Array[Array[File]] analysis_protocol_outputs = CreateIntermediateOptimusAdapters.analysis_protocol_outputs
    Array[Array[File]] links_outputs = CreateIntermediateOptimusAdapters.links_outputs
    Array[Array[File]] loom_file_descriptor_outputs = CreateIntermediateOptimusAdapters.loom_file_descriptor_outputs
    Array[Array[File]?] bam_file_descriptor_outputs = CreateIntermediateOptimusAdapters.bam_file_descriptor_outputs
  }

}