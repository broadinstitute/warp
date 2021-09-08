version 1.0

import "../tasks/AdapterTasks.wdl" as Tasks

# Making this a subworkflow will (hopefully) allow us to get the pertinent task level iput for the metadata to include in the analysis_process file

workflow MergeOptimusLooms {
  meta {
    description: "Creates a merged project level matrix for Optimus"
    allowNestedInputs: true
  }

  input {
    Array[File] output_looms
    String library
    String species
    String organ
    String project_id
    String project_name
    String output_basename
  }

  String pipeline_version = "1.0.0"

  call Tasks.MergeLooms as MergeLooms {
    input:
      output_looms = output_looms,
      library = library,
      species = species,
      organ = organ,
      project_id = project_id,
      project_name = project_name,
      output_basename = output_basename
  }

  call Tasks.GetOptimusPipelineVersion {
    input:
      pipeline_version = pipeline_version
  }

  output {
    File project_loom = MergeLooms.project_loom
    String pipeline_version_string = GetOptimusPipelineVersion.pipeline_version_string
  }
}

