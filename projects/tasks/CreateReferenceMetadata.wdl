version 1.0

import "../tasks/AdapterTasks.wdl" as Tasks

workflow CreateReferenceMetadata {
  meta {
    description: "Creates file_descriptor and reference_file metadata objects"
    allowNestedInputs: true
  }

  input {
    Array[String] reference_fastas
    String pipeline_type
    String version_timestamp
    String input_type
    String species
  }

  call Tasks.CheckInput as CheckReferences {
    input:
      input_array = reference_fastas,
      input_type = "reference",
      illegal_characters = "; ="
  }

  call Tasks.GetReferenceDetails {
    input:
      ref_fasta = CheckReferences.output_string,
      species = species
  }

  call Tasks.GetCloudFileCreationDate {
    input:
      file_path = CheckReferences.output_string
  }

  call Tasks.GetFileDescriptor as CreateReferenceFileDescriptor {
    input:
      file_path = CheckReferences.output_string,
      file_path_string = CheckReferences.output_string,
      input_uuid = CheckReferences.output_string, # Reference files do not have a unique id, so the file path is hashed to create an id
      pipeline_type = pipeline_type,
      creation_time = GetCloudFileCreationDate.creation_date,
      version_timestamp = version_timestamp
  }

  call Tasks.GetReferenceFileMetadata {
    input:
      file_path = CheckReferences.output_string,
      input_uuid = CheckReferences.output_string, # Reference files do not have a unique id, so the file path is hashed to create an id
      genus_species = species,
      assembly_type = GetReferenceDetails.assembly_type,
      pipeline_type = pipeline_type,
      ncbi_taxon_id = GetReferenceDetails.ncbi_taxon_id,
      reference_type = GetReferenceDetails.reference_type,
      version_timestamp = version_timestamp,
      reference_version = GetReferenceDetails.reference_version
  }

  output {
    Array[File] reference_file_descriptor_outputs = CreateReferenceFileDescriptor.file_descriptor_outputs
    String reference_file_uuid = GetReferenceFileMetadata.reference_file_uuid
    Array[File] reference_metadata_outputs = GetReferenceFileMetadata.reference_metadata_outputs
    String reference_fasta = CheckReferences.output_string
  }
}

