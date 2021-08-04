version 1.0

import "../../projects/tasks/AdapterTasks.wdl" as Tasks

workflow CreateReferenceMetadata {
  meta {
    description: "Creates file_descriptor and reference_file metadata objects"
    allowNestedInputs: true
  }

  input {
    Array[File] reference_fastas
    String Species
    String pipeline_type
    String workflow_version

  }


  # Parse the refernce information form the cromwell metadata and confirm that all workflows used the same reference
  call Tasks.CheckReferences {
    input:
      reference_fastas = reference_fastas
  }

  File reference_fasta = reference_fasta[0]

  call Tasks.GetReferenceInfo {
    input:
      reference_file = reference_fasta,
      species = species
  }

  call Tasks.CreateFileDescriptor as CreateReferenceFileDescriptor {
    input:
      reference_file = GetRefernce.reference
  }


  call Tasks.GetReferenceFileMetadata {
    input:
      file_path = reference_fasta,
      input_uuid = reference_fasta,  # Same as the file_path ?!
      genus_species = species,
      assembly_type = GetReferenceInfo.assembly,
      pipeline_type = pipeline_type,
      ncbi_taxon_id = GetReferenceInfo.ncbi_taxon_id,
      reference_type = GetReferenceInfo.reeference_type,
      version_timestamp = workflow_version,
      reference_version = GetReferenceInfo.reference_version
  }


  output {
    File reference_file_descriptor = CreateReferenceFileDescriptor.file_descriptor
    File reference_file_metadata = CreateReferenceFileDescriptor.metadata
  }
}

