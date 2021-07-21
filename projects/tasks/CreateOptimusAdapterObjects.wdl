version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing

workflow CreateOptimusAdapterObjects {
  meta {
    description: "Creates json objects for indexing HCA 10x analysis data"
    allowNestedInputs: true
  }

  input {
    Array[File] bams
    Array[File] looms
    String library
    String species
    String

    String cromwell_url = "https://api.firecloud.org/"

  }

  # Crerate  the project level loom
  call MergeLooms {
      input:
        # Fill in input for subworkflow
    }

  # Create the adapter json objects for each intermediate output
  scatter(idx in range(length(looms))) {
    File loom = looms[idx]
    File bam = bams[idx]

    call GetMetadata {
      input:
        output_path = loom
        cromwell_url = cromwell_url
        include_subworkflows = false # TODO: do we need subworkflows???
        include_keys =
    }

    call GetAnalysisFileMetadata{
      input:
    }
  }

  # Create the adapter json objects for the project matrtix




  output {
    Array[File] metadata_json = GetMetadata.metadata
  }
}

