version 1.0

import "../optimus/CreateOptimusAdapterObjects.wdl" as CreateOptimusObjects
import "../optimus/MergeOptimusLooms.wdl" as MergeLooms
import "../tasks/AdapterTasks.wdl" as Tasks
import "../tasks/CreateReferenceMetadata.wdl" as CreateReferenceMetadata


workflow CreateOptimusAdapterMetadata {
  meta {
    description: "Creates json objects for indexing HCA analysis data"
    allowNestedInputs: true
  }

  input {
    Array[File] output_bams
    Array[File] output_looms
    Array[String] input_ids #sequencing_process_provenance_document_id
    Array[String] fastq_1_uuids #array of space separated strings
    Array[String] fastq_2_uuids #array of space separated strings
    Array[String]? fastq_i1_uuids #array of space separated strings

    # These values come in as arrays from Terra, but should be populated with a single value (which may be repeated)
    Array[String] all_libraries
    Array[String] all_organs
    Array[String] all_species
    Array[String] all_project_ids
    Array[String] all_project_names
    
    # Flag which tells data importers if this run is a re-run of a previously processed project
    Boolean is_update

    String output_basename

    String cromwell_url = "https://api.firecloud.org/"
    String staging_area = "gs://broad-dsp-monster-hca-prod-lantern/"
    String pipeline_type = "Optimus"
    String workspace_bucket # gs:// path associated with the terra workspace

    String? version_timestamp # default behavior is to retrieve this from workspace_bucket metadata
  }

  ########################## Set up Inputs ##########################
  # version of this pipeline
  String pipeline_version = "1.0.0"

  # Get the version timestamp which is the creation date of the staging bucket
  call Tasks.GetBucketCreationDate as GetVersionTimestamp {
    input:
      bucket_path = workspace_bucket
  }

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
      illegal_characters = "; ="
  }

  String library = CheckLibrary.output_string
  String organ = CheckOrgan.output_string
  String species = CheckSpecies.output_string
  String project_id = CheckProjectID.output_string
  String project_name = CheckProjectName.output_string

  # Default should be timestamp from bucket creation but can be overwritten if data is updated after import
  String workspace_version_timestamp = select_first([version_timestamp, GetVersionTimestamp.version_timestamp])

  # Build staging bucket
  # When validating the staging bucket it can't end with '/'
  String staging_bucket = staging_area + project_id + "/"
  String staging_bucket_validation = staging_area + project_id

  String project_stratum_string = "project=" + project_id + ";library=" + library + ";species=" + species + ";organ=" + organ

  if (false) {
    String none = "None"
  }

  ########################## Get Optimus Metadata Files ##########################
  scatter (idx in range(length(output_looms))) {
    String? fastq_i1_uuid = if defined(fastq_i1_uuids) then select_first([fastq_i1_uuids])[idx] else none
    call CreateOptimusObjects.CreateOptimusAdapterObjects as CreateIntermediateOptimusAdapters {
      input:
        bam = output_bams[idx],
        loom = output_looms[idx],
        input_id = input_ids[idx],
        process_input_ids = select_all([fastq_1_uuids[idx],fastq_2_uuids[idx], fastq_i1_uuid]),
        project_id = project_id,
        version_timestamp = workspace_version_timestamp,
        cromwell_url = cromwell_url,
        is_project_level = false
    }
  }

  # store variable resulting from intermediate run
  Array[File] intermediate_links = flatten(CreateIntermediateOptimusAdapters.links_outputs)
  Array[File] intermediate_analysis_process_objects = flatten(CreateIntermediateOptimusAdapters.analysis_process_outputs)
  Array[File] intermediate_analysis_protocol_objects = flatten(CreateIntermediateOptimusAdapters.analysis_protocol_outputs)
  Array[File] intermediate_analysis_file_objects = flatten(CreateIntermediateOptimusAdapters.analysis_file_outputs)
  Array[File] intermediate_loom_descriptor_objects = flatten(CreateIntermediateOptimusAdapters.loom_file_descriptor_outputs)
  Array[File] intermediate_bam_descriptor_objects = flatten(select_all(CreateIntermediateOptimusAdapters.bam_file_descriptor_outputs))

  call CreateReferenceMetadata.CreateReferenceMetadata {
    input:
      reference_fastas = CreateIntermediateOptimusAdapters.reference_fasta,
      species = species,
      pipeline_type = pipeline_type,
      version_timestamp = workspace_version_timestamp,
      input_type = "reference"
  }

  Array[File] reference_fasta_array = [CreateReferenceMetadata.reference_fasta]

  # Merge all intermediate run looms to a single project level loom
  call MergeLooms.MergeOptimusLooms {
    input:
      output_looms = output_looms,
      library = library,
      species = species,
      organ = organ,
      project_id = project_id,
      project_name = project_name,
      output_basename = output_basename
  }

  Array[File] project_loom_array = [MergeOptimusLooms.project_loom]

  # Get all of the intermediate loom file
  call Tasks.GetProjectLevelInputIds {
    input:
      intermediate_analysis_files = flatten(CreateIntermediateOptimusAdapters.analysis_file_outputs)
  }

  # Create the project level objects based on the intermediate looms and the final merged loom
  call CreateOptimusObjects.CreateOptimusAdapterObjects as CreateProjectOptimusAdapters {
    input:
      loom = MergeOptimusLooms.project_loom,
      process_input_ids = [GetProjectLevelInputIds.process_input_uuids],
      input_id = project_stratum_string,
      project_id = project_id,
      version_timestamp = workspace_version_timestamp,
      cromwell_url = cromwell_url,
      is_project_level = true,
      reference_file_fasta = CreateIntermediateOptimusAdapters.reference_fasta[0],
      pipeline_version = MergeOptimusLooms.pipeline_version_string
  }

  call Tasks.CreateStagingAreaFile as CreateStagingAreaFile {
    input:
      is_update = is_update
  }

  # store variable resulting from project run
  Array[File] project_links = CreateProjectOptimusAdapters.links_outputs
  Array[File] project_analysis_process_objects = CreateProjectOptimusAdapters.analysis_process_outputs
  Array[File] project_analysis_protocol_objects = CreateProjectOptimusAdapters.analysis_protocol_outputs
  Array[File] project_analysis_file_objects = CreateProjectOptimusAdapters.analysis_file_outputs
  Array[File] project_loom_descriptor_objects = CreateProjectOptimusAdapters.loom_file_descriptor_outputs

  ########################## Copy Files to Staging Bucket ##########################
  Array[File] links_objects = flatten([intermediate_links, project_links])
  Array[File] analysis_file_descriptor_objects = flatten([intermediate_loom_descriptor_objects, intermediate_bam_descriptor_objects, project_loom_descriptor_objects])
  Array[File] analysis_file_metadata_objects = flatten([intermediate_analysis_file_objects, project_analysis_file_objects])
  Array[File] analysis_process_objects = flatten([intermediate_analysis_process_objects, project_analysis_process_objects])
  Array[File] analysis_protocol_objects = flatten([intermediate_analysis_protocol_objects, project_analysis_protocol_objects])
  Array[File] reference_metadata_objects = CreateReferenceMetadata.reference_metadata_outputs
  Array[File] reference_file_descriptor_objects = CreateReferenceMetadata.reference_file_descriptor_outputs
  Array[File] data_objects = flatten([reference_fasta_array, project_loom_array, output_bams, output_looms])
  File is_update_file = CreateStagingAreaFile.is_update_file

  call Tasks.CopyToStagingBucket as CopyToStagingBucket {
    input:
      staging_bucket = staging_bucket,
      links_objects = links_objects,
      analysis_file_descriptor_objects = analysis_file_descriptor_objects,
      analysis_file_metadata_objects = analysis_file_metadata_objects,
      analysis_process_objects = analysis_process_objects,
      analysis_protocol_objects = analysis_protocol_objects,
      reference_metadata_objects = reference_metadata_objects,
      reference_file_descriptor_objects = reference_file_descriptor_objects,
      data_objects = data_objects,
      is_update_file = is_update_file
  }

  # Only validate the staging bucket after files have been copied otherwise it will fail
  # The 'done' flag is a hack to make this depend on the completion of CopyToStagingBucket
  call Tasks.ValidateStagingArea {
    input:
      staging_area = staging_bucket_validation,
      done = CopyToStagingBucket.done
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
    File output_is_update_file = is_update_file
  }
}

