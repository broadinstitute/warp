version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing

workflow OptimusPostProcessing {
  meta {
    description: "Creates a combined matrix and the json files corresponding to the combined matrix for the HCA MVP"
    allowNestedInputs: true
  }

  input {
    Array[File] library_looms
    Array[File] analysis_file_jsons
    Array[File] links_jsons
    Array[String] library
    Array[String] species
    Array[String] organ
    String project_id
    String project_name
    String output_basename
    String staging_area = "gs://broad-dsp-monster-hca-prod-lantern/"
    String version_timestamp = "2021-05-24T12:00:00.000000Z"
  }


  # version of this pipeline
  String pipeline_version = "1.0.0"

  String project_stratum_string = "project=" + project_id + ";library=" + library[0] + ";species=" + species[0] + ";organ=" + organ[0]

  # Build staging bucket
  String staging_bucket = staging_area + project_id + "/staging/"

  call PostProcessing.CheckMetadata {
      input:
        library = library,
        species = species,
        organ = organ
  }

  call PostProcessing.MergeLooms {
    input:
      library_looms = library_looms,
      library = library[0],
      species = species[0],
      organ = organ[0],
      project_id = project_id,
      project_name = project_name,
      output_basename = output_basename
  }

  call PostProcessing.GetInputMetadata {
    input:
      analysis_file_jsons = analysis_file_jsons,
      output_basename = output_basename
  }

  call PostProcessing.CreateAdapterJson {
    input:
      project_loom = MergeLooms.project_loom,
      project_id = project_id,
      input_metadata_json = GetInputMetadata.input_metadata_json,
      project_stratum_string = project_stratum_string,
      staging_bucket = staging_bucket,
      version_timestamp = version_timestamp,
      pipeline_version = "optimus_post_processing_v" + pipeline_version
  }

  output {
    File project_loom = MergeLooms.project_loom
    Array[File] json_adapter_files = CreateAdapterJson.json_adapter_files
  }
}

