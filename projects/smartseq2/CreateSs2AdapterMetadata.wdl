version 1.0

import "../smartseq2/CreateSs2AdapterObjects.wdl" as CreateSs2Objects
import "../hca_mvp/tasks/AdapterTasks.wdl" as Tasks
import "../hca_mvp/tasks/CreateReferenceMetadata.wdl" as CreateReferenceMetadata


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
    Array[String] fastq_2_uuids #array of space separated strings
    Array[String]? fastq_i1_uuids #array of space separated strings

    # These values come in as arrays from Terra, but should be populated with a single value (which may be repeated)
    Array[String] all_libraries
    Array[String] all_organs
    Array[String] all_species
    Array[String] all_project_ids
    Array[String] all_project_names

    String cromwell_url = "https://firecloud-orchestration.dsde-dev.broadinstitute.org"
    String staging_area = "gs://fc-b4648544-9363-4a04-aa37-e7031c078a67/"
    String version_timestamp
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

  call Tasks.GetSs2PipelineVersion as CheckPipelineVersion {
    input:
      pipeline_version = pipeline_version
  }

  String library = CheckLibrary.output_string
  String organ = CheckOrgan.output_string
  String species = CheckSpecies.output_string
  String project_id = CheckProjectID.output_string
  String project_name = CheckProjectName.output_string


  # Build staging bucket
  String staging_bucket = staging_area + project_id + "/staging/"
  String project_stratum_string = "project=" + project_id + ";library=" + library + ";species=" + species + ";organ=" + organ

  if (false) {
    String none = "None"
  }

  call Tasks.GetCromwellMetadata {
    input:
      output_path = output_loom,
      cromwell_url = cromwell_url,
      include_subworkflows = true
  }

  call Tasks.ParseCromwellMetadata {
    input:
      cromwell_metadata = GetCromwellMetadata.metadata,
      pipeline_type = "SS2"
  }

  ########################## Get Ss2 Metadata Files ##########################
  scatter (idx in range(length(output_bams))) {
    String? fastq_i1_uuid = if defined(fastq_i1_uuids) then select_first([fastq_i1_uuids])[idx] else none
    call CreateSs2Objects.CreateSs2AdapterObjects as CreateIntermediateSs2Adapters {
      input:
        bam = output_bams[idx],
        bai = output_bais[idx],
        loom = output_loom,
        input_id = input_ids[idx],
        process_input_ids = select_all([fastq_1_uuids[idx],fastq_2_uuids[idx], fastq_i1_uuid]),
        project_id = project_id,
        version_timestamp = version_timestamp,
        pipeline_version = "ss2",
        reference_file_fasta = ParseCromwellMetadata.ref_fasta,
        metadata = GetCromwellMetadata.metadata,
        is_project_level = false
    }
  }

  # store variable resulting from intermediate run
  # Array[File] intermediate_links = flatten(CreateIntermediateSs2Adapters.links_outputs)
  Array[File] intermediate_analysis_process_objects = flatten(CreateIntermediateSs2Adapters.analysis_process_outputs)
  Array[File] intermediate_analysis_protocol_objects = flatten(CreateIntermediateSs2Adapters.analysis_protocol_outputs)
  Array[File] intermediate_analysis_file_objects = flatten(CreateIntermediateSs2Adapters.analysis_file_outputs)
  Array[File] intermediate_loom_descriptor_objects = flatten(select_all(CreateIntermediateSs2Adapters.loom_file_descriptor_outputs))
  Array[File] intermediate_bam_descriptor_objects = flatten(select_all(CreateIntermediateSs2Adapters.bam_file_descriptor_outputs))
  Array[File] intermediate_bai_descriptor_objects = flatten(select_all(CreateIntermediateSs2Adapters.bai_file_descriptor_outputs))


  # should this pipeline_type be Smartseq2_Multisample or SS2? does it matter?
  call CreateReferenceMetadata.CreateReferenceMetadata {
    input:
      reference_fastas = [ParseCromwellMetadata.ref_fasta],
      species = species,
      pipeline_type = "SS2",
      version_timestamp = version_timestamp,
      input_type = "reference"
  }

  Array[File] reference_fasta_array = [CreateReferenceMetadata.reference_fasta]

  # Get all of the intermediate loom file
  call Tasks.GetProjectLevelInputIds {
    input:
      intermediate_analysis_files = flatten(CreateIntermediateSs2Adapters.analysis_file_outputs)
  }

  # Create the project level objects based on the intermediate looms and the final merged loom
  call CreateSs2Objects.CreateSs2AdapterObjects as CreateProjectSs2Adapters {
    input:
      loom = output_loom,
      process_input_ids = [GetProjectLevelInputIds.process_input_uuids],
      input_id = project_stratum_string,
      project_id = project_id,
      version_timestamp = version_timestamp,
      is_project_level = true,
      reference_file_fasta = ParseCromwellMetadata.ref_fasta,
      pipeline_version = "ss2",
      metadata = GetCromwellMetadata.metadata
  }

  # store variable resulting from project run
  # Array[File] project_links = CreateProjectSs2Objects.links_outputs // TODO create large links file
  Array[File] project_analysis_process_objects = CreateProjectSs2Adapters.analysis_process_outputs
  Array[File] project_analysis_protocol_objects = CreateProjectSs2Adapters.analysis_protocol_outputs
  Array[File] project_analysis_file_objects = CreateProjectSs2Adapters.analysis_file_outputs
  Array[File] project_loom_descriptor_objects = flatten(select_all([CreateProjectSs2Adapters.loom_file_descriptor_outputs]))

  ########################## Copy Files to Staging Bucket ##########################
  # Array[File] links_objects = flatten([intermediate_links, project_links]) # TODO create large links file
  Array[File] analysis_file_descriptor_objects = flatten([intermediate_loom_descriptor_objects, intermediate_bam_descriptor_objects, intermediate_bai_descriptor_objects, project_loom_descriptor_objects])
  Array[File] analysis_file_metadata_objects = flatten([intermediate_analysis_file_objects, project_analysis_file_objects])
  Array[File] analysis_process_objects = flatten([intermediate_analysis_process_objects, project_analysis_process_objects])
  Array[File] analysis_protocol_objects = flatten([intermediate_analysis_protocol_objects, project_analysis_protocol_objects])
  Array[File] reference_metadata_objects = CreateReferenceMetadata.reference_metadata_outputs
  Array[File] reference_file_descriptor_objects = CreateReferenceMetadata.reference_file_descriptor_outputs
  Array[File] data_objects = flatten([reference_fasta_array, [output_loom], output_bams, output_bais])

  call Tasks.CopyToStagingBucket {
    input:
      staging_bucket = staging_bucket,
      # links_objects = links_objects,
      analysis_file_descriptor_objects = analysis_file_descriptor_objects,
      analysis_file_metadata_objects = analysis_file_metadata_objects,
      analysis_process_objects = analysis_process_objects,
      analysis_protocol_objects = analysis_protocol_objects,
      reference_metadata_objects = reference_metadata_objects,
      reference_file_descriptor_objects = reference_file_descriptor_objects,
      data_objects = data_objects
  }

  output {
    # Array[File] output_links_objects = links_objects
    Array[File] output_analysis_file_descriptor_objects = analysis_file_descriptor_objects
    Array[File] output_analysis_file_metadata_objects = analysis_file_metadata_objects
    Array[File] output_analysis_process_objects = analysis_process_objects
    Array[File] output_analysis_protocol_objects = analysis_protocol_objects
    Array[File] output_reference_metadata_objects = reference_metadata_objects
    Array[File] output_reference_file_descriptor_objects = reference_file_descriptor_objects
    Array[File] output_data_objects = data_objects
  }
}

