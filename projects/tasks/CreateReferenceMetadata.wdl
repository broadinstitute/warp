version 1.0

import "../../projects/tasks/AdapterTasks.wdl" as Tasks

workflow CreateReferenceMetadata {
  meta {
    description: "Creates file_descriptor and reference_file metadata objects"
    allowNestedInputs: true
  }

  input {
    Array[File] optimus_metadata_jsons
  }


  # Parse the refernce information form the cromwell metadata and confirm that all workflows used the same reference
  call Tasks.GetRefernce {
    input:
      optimus_metadata_jsons = optimus_metadata_jsons
  }


  call Tasks.CreateFileDescriptor as CreateReferenceFileDescriptor {
    input:
      reference_file = GetRefernce.reference
  }


  call Tasks.CreateReferenceFileMetadata {
    input:
      reference_file = GetRefernce.reference
  }


  output {
    File reference_file_descriptor = CreateReferenceFileDescriptor.file_descriptor
    File reference_file_metadata = CreateReferenceFileDescriptor.metadata
  }
}

