version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing

workflow CreateAdapterMetadata {
  meta {
    description: "Creates json objects for indexing HCA analysis data"
    allowNestedInputs: true
  }

  input {
    Array[File] output_bams
    Array[File] output_looms
    Array[File]? output_bais
    Array[String] input_ids

    # These values come in as arrays froom Terra, but should be populated with a single value (which may be repeated)
    Array[String] all_libraries
    Array[String] all_species
    Array[String] all_organs
    Array[String] all_project_ids
    Array[String] all_project_names

    String output_basename
    String staging_area = "gs://broad-dsp-monster-hca-prod-lantern/"
    String version_timestamp = "2021-05-24T12:00:00.000000Z"
  }

  # version of this pipeline
  String pipeline_version = "1.0.0"

  Boolean is_SS2 = false # TODO: check an input value to determine if this is SS2
  Boolean is_Optimus = true # TODO: check an input value to determine if this is Optimus (leaving this flexible for additional data types if needed

  call CheckStratumMetadata {
    input:
      library = all_libraries,
      species = all_species,
      organ = all_organs,
      project_id = all_project_ids,
      project_name = all_project_names
  }

  String library = all_libraries[0]
  String species = all_species[0]
  String organ = all_organs[0]
  String project_id = all_project_ids[0]
  String project_name = all_project_names[0]

  # Build staging bucket
  String staging_bucket = staging_area + project_id + "/staging/"
  String project_stratum_string = "project=" + project_id + ";library=" + library + ";species=" + species + ";organ=" + organ

  # Split out subworkflows based on datatype
  if (is_Optimus) {

      call CreateOptimusAdapterObjects as GetAdapterObjects{
        input:
          bam = output_bams[idx],
          loom = output_looms[idx],
          input_id = input_ids
    }
  }

  if (is_SS2) {
    call CreateSS2AdapterObjects as GetAdapterObjects{
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

