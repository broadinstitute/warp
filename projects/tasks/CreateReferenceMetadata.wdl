version 1.0

import "../../projects/tasks/AdapterTasks.wdl" as Tasks

workflow CreateReferenceMetadata {
  meta {
    description: "Creates file_descriptor and reference_file metadata objects"
    allowNestedInputs: true
  }

  input {
    Array[String] reference_fastas
    String pipeline_type
    String version_timestamp

    String species
    String pipeline_type
    String workflow_version

  }

  call Tasks.CheckReferences {
    input:
      reference_fastas = reference_fastas
  }

  String reference_fasta = reference_fastas[0]

  call Tasks.GetReferenceDetails {
    input:
      reference_file = reference_fasta,
      species = species
  }

  call Tasks.GetCloudFileCreationDate {
    input:
      file_path = reference_fasta
  }

  call Tasks.CreateFileDescriptor as CreateReferenceFileDescriptor {
    input:
      file_path = reference_fasta,
      file_path_string = reference_fasta,
      input_uuid = reference_fasta, # Reference files do not have a unique id, so the file path is hashed to create an id
      pipeline_type = pipeline_type,
      creation_time = GetCloudFileCreationDate.creation_date,
      version_timestamp = version_timestamp
  }

  call Tasks.GetReferenceFileMetadata {
    input:
      file_path = reference_fasta,
      input_uuid = reference_fasta, # Reference files do not have a unique id, so the file path is hashed to create an id
      genus_species = species,
      assembly_type = GetReferenceDetails.assembly_type,
      pipeline_type = pipeline_type,
      ncbi_taxon_id = GetReferenceDetails.ncbi_taxon_id,
      reference_type = GetReferenceDetails.reference_type,
      version_timestamp = version_timestamp,
      reference_version = GetReferenceInfo.reference_version # TODO: what is this supposed to be??
  }


  output {
    File reference_file_descriptor = CreateReferenceFileDescriptor.file_descriptor
    File reference_file_metadata = CreateReferenceFileDescriptor.metadata
  }
}

