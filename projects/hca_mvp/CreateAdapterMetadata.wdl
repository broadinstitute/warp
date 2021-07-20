version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing

workflow CreateAdapterMetadata {
  meta {
    description: "Creates json objects for indexing HCA analysis data"
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

  Boolean is_SS2 = false # TODO: check an input value to determine if this is SS2
  Boolean is_Optimus = true # TODO: check an input value to determine if this is Optimus (leaving this flexible for additional data types if needed)

  # version of this pipeline
  String pipeline_version = "1.0.0"

  String project_stratum_string = "project=" + project_id + ";library=" + library[0] + ";species=" + species[0] + ";organ=" + organ[0]

  # Build staging bucket
  String staging_bucket = staging_area + project_id + "/staging/"

  # Split out subworkflows based on datatype
  if (is_SS2) {
    call CreateSS2AdapterObjects as GetAdapterObjects{
      input:
        # Fill in input for subworkflow
    }
  }

  if (is_Optimus) {
    call MergeLooms {
      input:
        # Fill in input for subworkflow
    }
    call CreateOptimusAdapterObjects as GetAdapterObjects{
      input:
        # Fill in input for subworkflow
    }
  }


  call CopyToStagingBucket {
    input:
      Array[File] links_objects = GetAdapterObjects.json_adapters,
      Array[File] desscriptor_objects = ,
      Array[File]
      Array[File] data = data

  }


  output {
    Array[File] json_adapters = GetAdapterObjects.json_adapters
  }
}

