version 1.0

import "../../projects/tasks/CreateOptimusAdapterObjects.wdl" as CreateOptimusObjects
import "../../projects/tasks/MergeOptimusLooms.wdl" as MergeLooms
import "../../projects/tasks/AdapterTasks.wdl" as Tasks
import "../../projects/tasks/CreateReferenceMetadata.wdl" as CreateReferenceMetadata


workflow CreateAdapterMetadata {
  meta {
    description: "Creates json objects for indexing HCA analysis data"
    allowNestedInputs: true
  }

  input {
    Array[File] output_bams
    Array[File] output_looms
    Array[File]? output_bais
    Array[String] input_ids #sequencing_process_provenance_document_id
    Array[String] fastq_1_uuids
    Array[String] fastq_2_uuids
    Array[String] fastq_i1_uuids = []

    # These values come in as arrays from Terra, but should be populated with a single value (which may be repeated)
    Array[String] all_libraries
    Array[String] all_species
    Array[String] all_organs
    Array[String] all_project_ids
    Array[String] all_project_names

    String output_basename
    String cromwell_url = "https://api.firecloud.org/"
    String staging_area = "gs://broad-dsp-monster-hca-prod-lantern/"
    String version_timestamp = "2021-05-24T12:00:00.000000Z" # TODO should we hard code this?
  }

  ########################## Set up Inputs ##########################
  # version of this pipeline
  String pipeline_version = "1.0.0"

  Boolean is_SS2 = false # TODO: check an input value to determine if this is SS2
  Boolean is_Optimus = true # TODO: check an input value to determine if this is Optimus (leaving this flexible for additional data types if needed

  call Tasks.CheckStratumMetadata {
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

  Array[String] fastq_uuids = flatten([fastq_1_uuids, fastq_2_uuids, fastq_i1_uuids])

  ########################## Get Optimus Metadata Files ##########################
  if (is_Optimus) {
    scatter (idx in range(length(looms))) {
      call CreateOptimusObjects.CreateOptimusAdapterObjects as GetIntermediateOptimusAdapters {
        input:
          bam = bams[idx],
          loom = looms[idx],
          input_id = input_ids[idx],
          library = library,
          species = species,
          organ = organ,
          project_id = project_id,
          project_name = project_name,
          version_timestamp = version_timestamp,
          cromwell_url = cromwell_url,
          project_level = false

      }
    }
    call CreateReferenceMetadata.CreateReferenceMetadata as CreateReferenceMetadata {
      input:
        reference_fastas = reference_fastas,
        species = species,
        pipeline_type = pipeline_type,
        workflow_version = workflow_version
    }
    call MergeLooms.MergeOptimusLooms as MergeLooms {
      input:
        output_looms = output_looms,
        library = library,
        species = species,
        organ = organ,
        project_id = project_id,
        project_name = project_name,
        output_basename = output_basename
    }
    # get adapters for merged matrix
    call CreateOptimusObjects.CreateOptimusAdapterObjects as GetProjectOpitmusAdapters {
      input:
        loom = MergeLooms.project_loom,
        input_id = project_stratum_string,
        library = library,
        species = species,
        organ = organ,
        project_id = project_id,
        project_name = project_name,
        project_stratum_string = project_stratum_string
      }
    }
  }

  ########################## Get SS2 Metadata Files ###########################
  if (is_SS2) {
    call CreateSS2Objects as GetAdapterObjects{
      input:
        # Fill in input for subworkflow
    }
  }

  ########################## Copy Files to Staging Bucket ##########################
  call CopyToStagingBucket {
    input:
      Array[File] links_objects = GetAdapterObjects.json_adapters,
      Array[File] desscriptor_objects = ,
      Array[File]
      Array[File] data = data
      String? cache_invalidate

  }


  output {
    Array[File] json_adapters = GetAdapterObjects.json_adapters
  }
}

