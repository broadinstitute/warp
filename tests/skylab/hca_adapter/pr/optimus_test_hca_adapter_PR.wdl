version 1.0

import "../../../../projects/optimus/CreateOptimusAdapterMetadata.wdl" as test_target_adapter
import "../../../../tests/skylab/hca_adapter/pr/ValidateHcaAdapter.wdl" as checker_adapter

# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestOptimusHcaAdapter {
  input {
    # hca optimus inputs 
    Array[File] output_bams
    Array[File] output_looms
    Array[String] input_ids
    Array[String] fastq_1_uuids
    Array[String] fastq_2_uuids
    Array[String]? fastq_i1_uuids

    Array[String] all_libraries
    Array[String] all_organs
    Array[String] all_species
    Array[String] all_project_ids
    Array[String] all_project_names

    String output_basename
    String cromwell_url = "https://api.firecloud.org/"
    String staging_area = "gs://broad-dsp-monster-hca-prod-lantern/"
    String pipeline_type = "Optimus"

    #optimus truth inputs
    File optimus_descriptors_analysis_file_intermediate_loom_json
    File optimus_metadata_analysis_file_intermediate_bam_json
    File optimus_metadata_reference_file_intermediate_json
    File optimus_metadata_analysis_protocol_file_intermediate_json
    File optimus_links_intermediate_loom_json
    File optimus_metadata_analysis_process_file_intermediate_json
    File optimus_descriptors_analysis_file_intermediate_bam_json
    File optimus_metadata_analysis_file_intermediate_loom_json
    File optimus_descriptors_analysis_file_intermediate_reference_json
    File optimus_descriptors_analysis_file_project_loom_json
    File optimus_links_project_loom_json
    File optimus_metadata_analysis_file_project_loom_json
    File optimus_metadata_analysis_process_project_loom_json
    File optimus_metadata_analysis_protocol_project_loom_json
  }


  call test_target_adapter.CreateOptimusAdapterMetadata as target_adapter {
    input:
      output_bams = output_bams,
      output_looms = output_looms,
      input_ids = input_ids,
      fastq_1_uuids = fastq_1_uuids,
      fastq_2_uuids = fastq_2_uuids,
      fastq_i1_uuids = fastq_i1_uuids,
      all_libraries = all_libraries,
      all_organs = all_organs,
      all_species = all_species,
      all_project_ids = all_project_ids,
      all_project_names = all_project_names,
      output_basename = output_basename,
      cromwell_url = cromwell_url,
      staging_area = staging_area,
      pipeline_type = pipeline_type
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_loom {
    input:
      test_json = target_adapter.output_analysis_file_descriptor_objects[0],
      truth_json = optimus_descriptors_analysis_file_intermediate_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_bam {
    input:
      test_json = target_adapter.output_analysis_file_descriptor_objects[1],
      truth_json = optimus_descriptors_analysis_file_intermediate_bam_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_reference {
    input:
      test_json = target_adapter.output_reference_file_descriptor_objects[0],
      truth_json = optimus_descriptors_analysis_file_intermediate_reference_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_links {
    input:
      test_json = target_adapter.output_links_objects[0],
      truth_json = optimus_links_intermediate_loom_json,
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_files_loom {
    input:
      test_json = target_adapter.output_analysis_file_metadata_objects[0],
      truth_json = optimus_metadata_analysis_file_intermediate_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_files_bam {
    input:
      test_json = target_adapter.output_analysis_file_metadata_objects[1],
      truth_json = optimus_metadata_analysis_file_intermediate_bam_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_process {
    input:
      test_json = target_adapter.output_analysis_process_objects[0],
      truth_json = optimus_metadata_analysis_process_file_intermediate_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_reference_file {
      input:
        test_json = target_adapter.output_reference_metadata_objects[0],
        truth_json = optimus_metadata_reference_file_intermediate_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_protocol {
    input:
      test_json = target_adapter.output_analysis_protocol_objects[0],
      truth_json = optimus_metadata_analysis_protocol_file_intermediate_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_descriptors {
    input:
      test_json = target_adapter.output_analysis_file_descriptor_objects[2],
      truth_json = optimus_descriptors_analysis_file_project_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_links {
    input:
      test_json = target_adapter.output_links_objects[1],
      truth_json = optimus_links_project_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_metadata {
    input:
      test_json = target_adapter.output_analysis_file_metadata_objects[2],
      truth_json = optimus_metadata_analysis_file_project_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_metadata_analysis_process {
    input:
      test_json = target_adapter.output_analysis_process_objects[1],
      truth_json = optimus_metadata_analysis_process_project_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_metadata_analysis_protocol {
    input:
      test_json = target_adapter.output_analysis_protocol_objects[1],
      truth_json = optimus_metadata_analysis_protocol_project_loom_json
  }

    output {
      Array[File] analysis_file = target_adapter.output_analysis_file_metadata_objects
      Array[File] analysis_process = target_adapter.output_analysis_process_objects
      Array[File] analysis_protocol = target_adapter.output_analysis_protocol_objects
      Array[File] analysis_output = target_adapter.output_data_objects
      File reference_genome = target_adapter.output_data_objects[0]
      Array[File] reference_genome_reference_file = target_adapter.output_reference_metadata_objects
      Array[File] reference_genome_descriptor = target_adapter.output_reference_file_descriptor_objects
      Array[File] analysis_file_descriptor = target_adapter.output_analysis_file_descriptor_objects
      Array[File] links = target_adapter.output_links_objects
    }
}