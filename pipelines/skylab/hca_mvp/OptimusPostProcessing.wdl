version 1.0

import "../../../tasks/skylab/PostProcessingOptimus.wdl" as PostProcessing

workflow OptimusPostProcessing {
  meta {
    description: "Adds additional metadata to each looom ouptut and create json file anf a project level matrix"
  }

  input {
    Array[File] library_looms
    Array[String] library
    Array[String] species
    Array[String] stage
    Array[String] organ
    String output_basename
  }

  # version of this pipeline

  String pipeline_version = "1.0.0"



  parameter_meta {
  }
  call PostProcessing.CheckMetadata {
      input:
        library = library,
        species = species,
        stage = stage,
        organ = organ
  }

  #call PostProcessing.CreateMetadataTsv {
  #  input:
  #    original_looms = original_looms,
  #    library = library,
  #    species = species,
  #    stage = stage,
  #    organ = organ,
  #    output_basename = output_basename
  #}
  #scatter (loom in original_looms) {
  #  call PostProcessing.IndividualLoomProcessing {
  #      input:
  #        original_library_loom = loom,
  #        library = CheckMetadata.library,
  #        species = CheckMetadata.species,
  #        stage = CheckMetadata.stage,
  #        organ = CheckMetadata.organ
  #    }
  #}

  #Array[File] library_looms = IndividualLoomProcessing.library_loom
  #Array[File] library_jsons = IndividualLoomProcessing.library_json

  call PostProcessing.MergeLooms {
    input:
      library_looms = library_looms,
      library = CheckMetadata.library,
      species = CheckMetadata.species,
      stage = CheckMetadata.stage,
      organ = CheckMetadata.organ,
      output_basename = output_basename
  }

  output {
    #Array[File] library_looms = library_looms
    #Array[File] library_jsons = library_jsons

    #File project_tsv = CreateMetadataTsv.metadata_tsv

    File project_loom = MergeLooms.merged_loom
    File project_json = MergeLooms.merged_json
  }

  meta {
      allowNestedInputs: true
  }
}

