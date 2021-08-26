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
    Array[String]? fastq_i1_uuids

    # These values come in as arrays from Terra, but should be populated with a single value (which may be repeated)
    Array[String] all_libraries
    Array[String] all_organs
    Array[String] all_species
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

  # Check inputs for multiple values or illegal characters
  call Tasks.CheckInput as CheckLibrary {
    input:
      input_array = all_libraries,
      input_type = "library",
      illegal_characters = "; ="
  }

  call Tasks.CheckInput as CheckOrgan {
    input:
      input_array = all_organs,
      input_type = "organ",
      illegal_characters = "; ="
  }

  call Tasks.CheckInput as CheckSpecies {
    input:
      input_array = all_species,
      input_type = "species",
      illegal_characters = "; ="
  }

  call Tasks.CheckInput as CheckProjectID {
    input:
      input_array = all_project_ids,
      input_type = "project_id",
      illegal_characters = "; ="
  }

    call Tasks.CheckInput as CheckProjectName {
    input:
      input_array = all_project_names,
      input_type = "project_name",
      illegal_characters = "; =" # should we include % in this list? # ultimately we should switch to a whitelist
  }

  String library = CheckLibrary.output_string
  String organ = CheckOrgan.output_string
  String species = CheckSpecies.output_string
  String project_id = CheckProjectID.output_string
  String project_name = CheckProjectName.output_string

  call Tasks.GetPipelineType {
    input:
      library = library
  }

  # Build staging bucket
  String staging_bucket = staging_area + project_id + "/staging/"
  String project_stratum_string = "project=" + project_id + ";library=" + library + ";species=" + species + ";organ=" + organ

  if (false) {
    String none = "None"
  }
  ########################## Get Optimus Metadata Files ##########################
  if (GetPipelineType.output_string == "Optimus") {
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
    call CreateReferenceMetadata.CreateReferenceMetadata as CreateReferenceMetadata {
      input:
        reference_fastas = CreateIntermediateOptimusAdapters.reference_fasta,
        species = species,
        pipeline_type = 'Optimus',
        version_timestamp = version_timestamp,
        input_type = "reference"
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


    call Tasks.GetProjectLevelInputIds {
      input:
        intermediate_analysis_files = flatten(CreateIntermediateOptimusAdapters.analysis_file_outputs)
    }

    # TODO: change the way the links file uses the project level input ids (need to accept strings instead of jsons)

    call CreateOptimusObjects.CreateOptimusAdapterObjects as CreateProjectOptimusAdapters {
      input:
        loom = MergeLooms.project_loom,
        process_input_ids = [GetProjectLevelInputIds.process_input_uuids],
        input_id = project_stratum_string,
        library = library,
        species = species,
        organ = organ,
        project_id = project_id,
        project_name = project_name,
        project_stratum_string = project_stratum_string,
        version_timestamp = version_timestamp,
        cromwell_url = cromwell_url,
        is_project_level = true,
        reference_file_fasta = CreateIntermediateOptimusAdapters.reference_fasta[0],
        pipeline_version = CreateIntermediateOptimusAdapters.pipeline_version_string[0]
    }
  }

  ########################## Get SS2 Metadata Files ###########################
  #if (GetPipelineType.output_string == "SS2") {
  #  call CreateSS2Objects as GetAdapterObjects{
  #    input:
  #      # Fill in input for subworkflow
  #  }
  #}

  ########################## Copy Files to Staging Bucket ##########################
    Array[File] links_objects = flatten(select_all([CreateIntermediateOptimusAdapters.links_outputs, CreateProjectOptimusAdapters.links_outputs]))
    Array[File] analysis_file_descriptor_objects = flatten(select_all([select_all([CreateIntermediateOptimusAdapters.loom_file_descriptor_outputs, CreateIntermediateOptimusAdapters.bam_file_descriptor_outputs]), CreateProjectOptimusAdapters.loom_file_descriptor_outputs]))
    #Array[File] analysis_file_descriptor_objects = flatten([CreateIntermediateOptimusAdapters.loom_file_descriptor_outputs, select_all([CreateIntermediateOptimusAdapters.bam_file_descriptor_outputs]), CreateProjectOptimusAdapters.loom_file_descriptor_outputs])
    Array[File] analysis_file_metadata_objects = flatten(select_all([CreateIntermediateOptimusAdapters.analysis_file_outputs, CreateProjectOptimusAdapters.analysis_file_outputs]))
    Array[File] analysis_process_objects = flatten(select_all([CreateIntermediateOptimusAdapters.analysis_process_outputs, CreateProjectOptimusAdapters.analysis_process_outputs]))
    Array[File] analysis_protocol_objects = flatten(select_all([CreateIntermediateOptimusAdapters.analysis_protocol_outputs, CreateProjectOptimusAdapters.analysis_protocol_outputs]))
    Array[File] reference_metadata_objects = select_first([CreateReferenceMetadata.reference_metadata_outputs])
    Array[File] reference_file_descriptor_objects = select_first([CreateReferenceMetadata.reference_file_descriptor_outputs])
    File? reference_fasta_file = CreateReferenceMetadata.reference_fasta
    Array[File] data_objects = flatten([select_all([output_bams, output_looms, reference_fasta_file, MergeLooms.project_loom])])

    call Tasks.CopyToStagingBucket {
      input:
        staging_bucket = staging_bucket,
        links_objects = links_objects,
        analysis_file_descriptor_objects = analysis_file_descriptor_objects,
        analysis_file_metadata_objects = analysis_file_metadata_objects,
        analysis_process_objects = analysis_process_objects,
        analysis_protocol_objects = analysis_protocol_objects,
        reference_metadata_objects = reference_metadata_objects,
        reference_file_descriptor_objects = reference_file_descriptor_objects,
        data_objects = data_objects
      }


  output {
    Array[File] output_links_objects = links_objects
    Array[File] output_analysis_file_descriptor_objects = analysis_file_descriptor_objects
    Array[File] output_analysis_file_metadata_objects = analysis_file_metadata_objects
    Array[File] output_analysis_process_objects = analysis_process_objects
    Array[File] output_analysis_protocol_objects = analysis_protocol_objects
    Array[File] output_reference_metadata_objects = reference_metadata_objects
    Array[File] output_reference_file_descriptor_objects = reference_file_descriptor_objects
    Array[File] output_data_objects = data_objects
  }
}

