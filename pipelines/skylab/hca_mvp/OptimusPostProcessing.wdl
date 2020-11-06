version 1.0

import "../../../tasks/skylab/OptimusPostProcessingTasks.wdl" as PostProcessing

workflow OptimusPostProcessing {
  meta {
    description: "Adds additional metadata to each looom ouptut and create json file anf a project level matrix"
  }

  input {
    Array[File] library_looms
    Array[String] library
    Array[String] species
    Array[String] organ
    String output_basename
  }

  # version of this pipeline
  String pipeline_version = "1.0.0"

  call PostProcessing.CheckMetadata {
      input:
        library = library,
        species = species,
        organ = organ
  }

  call PostProcessing.MergeLooms {
    input:
      library_looms = library_looms,
      library = library[0],
      species = species[0],
      organ = organ[0],
      output_basename = output_basename
  }

  call PostProcessing.CreateAdapterJson {
    input:
      project_loom = MergeLooms.project_loom,
      output_basename = output_basename
  }

  output {
    File project_loom = MergeLooms.project_loom
    File project_json = CreateAdapterJson.project_json
  }

  meta {
      allowNestedInputs: true
  }
}

