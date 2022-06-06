version 1.0

import "../smartseq2/CreateSs2AdapterObjects.wdl" as CreateSs2Objects
import "../tasks/AdapterTasks.wdl" as Tasks
import "../tasks/CreateReferenceMetadata.wdl" as CreateReferenceMetadata


workflow CreateSs2AdapterMetadata {
  meta {
    description: "Creates json objects for indexing HCA smartseq2 analysis data"
    allowNestedInputs: true
  }

  input {
    Array[File] output_bams
    Array[File] output_bais
    File output_loom
    Array[String] input_ids #sequencing_process_provenance_document_id
    Array[String] fastq_1_uuids #array of space separated strings
    Array[String]? fastq_2_uuids #array of space separated strings

    # These values come in as arrays from Terra, but should be populated with a single value (which may be repeated)
    Array[String] all_libraries
    Array[String] all_organs
    Array[String] all_species
    Array[String] all_project_ids
    Array[String] all_project_names

    # Flag which tells data importers if this run is a re-run of a previously processed project
    Boolean is_update

    String cromwell_url = "https://api.firecloud.org/"
    String staging_area = "gs://broad-dsp-monster-hca-prod-lantern/"
    String pipeline_type = "SS2"
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
  # When validating the staging bucket it can't end with '/
  String staging_bucket = staging_area + project_id + "/"
  String staging_bucket_validation = staging_area + project_id

  String project_stratum_string = "project=" + project_id + ";library=" + library + ";species=" + species + ";organ=" + organ

  call Tasks.GetCromwellMetadata {
    input:
      output_path = output_loom,
      cromwell_url = cromwell_url,
      include_subworkflows = true
  }

  call Tasks.ParseCromwellMetadata {
    input:
      cromwell_metadata = GetCromwellMetadata.metadata,
      pipeline_type = pipeline_type
  }

  # store variables from ParseCromwellMetadata
  String reference_fasta = ParseCromwellMetadata.ref_fasta
  String multi_sample_pipeline_version = ParseCromwellMetadata.pipeline_version
  String single_sample_pipeline_version = ParseCromwellMetadata.single_sample_pipeline_version

  ########################## Get Ss2 Metadata Files ##########################
  scatter (idx in range(length(output_bams))) {
    call CreateSs2Objects.CreateSs2AdapterObjects as CreateIntermediateSs2Adapters {
      input:
        bam = output_bams[idx],
        bai = output_bais[idx],
        input_id = input_ids[idx],
        ss2_index = idx,
        version_timestamp = workspace_version_timestamp,
        pipeline_version = single_sample_pipeline_version,
        pipeline_type = pipeline_type,
        reference_file_fasta = reference_fasta,
        metadata = GetCromwellMetadata.metadata,
        is_project_level = false
    }
  }

  # store variable resulting from intermediate run
  Array[File] intermediate_analysis_process_objects = flatten(CreateIntermediateSs2Adapters.analysis_process_outputs)
  Array[File] intermediate_analysis_protocol_objects = flatten(CreateIntermediateSs2Adapters.analysis_protocol_outputs)
  Array[File] intermediate_analysis_file_objects = flatten(CreateIntermediateSs2Adapters.analysis_file_outputs)
  Array[File] intermediate_loom_descriptor_objects = flatten(select_all(CreateIntermediateSs2Adapters.loom_file_descriptor_outputs))
  Array[File] intermediate_bam_descriptor_objects = flatten(select_all(CreateIntermediateSs2Adapters.bam_file_descriptor_outputs))
  Array[File] intermediate_bai_descriptor_objects = flatten(select_all(CreateIntermediateSs2Adapters.bai_file_descriptor_outputs))


  call CreateReferenceMetadata.CreateReferenceMetadata {
    input:
      reference_fastas = [reference_fasta],
      species = species,
      pipeline_type = pipeline_type,
      version_timestamp = workspace_version_timestamp,
      input_type = "reference"
  }

  Array[File] reference_fasta_array = [CreateReferenceMetadata.reference_fasta]

  # Create the project level adapter objects based on project loom
  call CreateSs2Objects.CreateSs2AdapterObjects as CreateProjectSs2Adapters {
    input:
      loom = output_loom,
      input_id = project_stratum_string,
      version_timestamp = workspace_version_timestamp,
      is_project_level = true,
      reference_file_fasta = reference_fasta,
      pipeline_version = multi_sample_pipeline_version,
      pipeline_type = pipeline_type,
      metadata = GetCromwellMetadata.metadata
  }

  # store variable resulting from project run
  Array[File] project_analysis_process_objects = CreateProjectSs2Adapters.analysis_process_outputs
  Array[File] project_analysis_protocol_objects = CreateProjectSs2Adapters.analysis_protocol_outputs
  File analysis_file_outputs_json = CreateProjectSs2Adapters.analysis_file_outputs_json
  Array[File] project_analysis_file_objects = CreateProjectSs2Adapters.analysis_file_outputs
  Array[File] project_loom_descriptor_objects = flatten(select_all([CreateProjectSs2Adapters.loom_file_descriptor_outputs]))


  call Tasks.GetLinksFileMetadata {
    input:
      project_id = project_id,
      output_file_path = analysis_file_outputs_json,
      version_timestamp = workspace_version_timestamp,
      process_input_ids = input_ids,
      analysis_process_path = project_analysis_process_objects,
      analysis_protocol_path = project_analysis_protocol_objects,
      analysis_process_path_list = intermediate_analysis_process_objects,
      analysis_protocol_path_list = intermediate_analysis_protocol_objects,
      bam_array = output_bams,
      bai_array = output_bais,
      fastq1_array = fastq_1_uuids,
      fastq2_array = select_first([ fastq_2_uuids, [] ]),
      project_level = true,
      file_name_string = project_stratum_string,
      pipeline_type = pipeline_type
  }

  call Tasks.CreateStagingAreaFile as CreateStagingAreaFile {
    input:
      is_update = is_update
  }

  Array[File] project_links = GetLinksFileMetadata.links_outputs

  ########################## Copy Files to Staging Bucket ##########################
  Array[File] links_objects = project_links
  Array[File] analysis_file_descriptor_objects = flatten([intermediate_loom_descriptor_objects, intermediate_bam_descriptor_objects, intermediate_bai_descriptor_objects, project_loom_descriptor_objects])
  Array[File] analysis_file_metadata_objects = flatten([intermediate_analysis_file_objects, project_analysis_file_objects])
  Array[File] analysis_process_objects = flatten([intermediate_analysis_process_objects, project_analysis_process_objects])
  Array[File] analysis_protocol_objects = flatten([intermediate_analysis_protocol_objects, project_analysis_protocol_objects])
  Array[File] reference_metadata_objects = CreateReferenceMetadata.reference_metadata_outputs
  Array[File] reference_file_descriptor_objects = CreateReferenceMetadata.reference_file_descriptor_outputs
  Array[File] data_objects = flatten([reference_fasta_array, [output_loom], output_bams, output_bais])
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
