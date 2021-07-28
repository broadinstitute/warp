version 1.0

import "../../projects/tasks/OptimusPostProcessingTasks.wdl" as PostProcessing

workflow CreateOptimusAdapterObjects {
  meta {
    description: "Creates json objects for indexing HCA 10x analysis data"
    allowNestedInputs: true
  }

  input {
    File input_file
    File metadata
    String library
    String species
    String organ
    String project_id
    String project_name
    String project_stratum_string

    String cromwell_url = "https://api.firecloud.org/"

  }




  call GetAnalysisFileMetadata{
    input:
  }

  # Create the adapter json objects for the project matrtix




  output {
    Array[File] metadata_json = GetMetadata.metadata
  }
}

