version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing

workflow CreateOptimusAdapterObjects {
  meta {
    description: "Creates json objects for indexing HCA 10x analysis data"
    allowNestedInputs: true
  }

  input {

  }

  call GetMetadata {
    input:
  }

  call GetAnalysisFileMetadata{
      input:
  }


  output {
    Array[File] metadata_json = GetMetadata.metadata
  }
}

