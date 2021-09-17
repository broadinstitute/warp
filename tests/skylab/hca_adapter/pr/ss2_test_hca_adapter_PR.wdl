version 1.0

import "../../../../projects/smartseq2/CreateSs2AdapterMetadata.wdl" as test_target_adapter
import "../../../../tests/skylab/hca_adapter/pr/ValidateHcaAdapter.wdl" as checker_adapter

# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestSs2HcaAdapter {
  input {
    # hca ss2 inputs 
    Array[File] output_bams
    Array[File] output_bais
    File output_loom
    Array[String] input_ids
    Array[String] fastq_1_uuids
    Array[String] fastq_2_uuids
    Array[String]? fastq_i1_uuids

    Array[String] all_libraries
    Array[String] all_organs
    Array[String] all_species
    Array[String] all_project_ids
    Array[String] all_project_names

    String workspace_bucket

    # ss2 truth inputs
    File ss2_metadata_analysis_file_intermediate_bam_json
    File ss2_metadata_analysis_file_intermediate_bai_json
    File ss2_metadata_analysis_file_project_loom_json
    File ss2_descriptors_analysis_file_intermediate_reference_json
    File ss2_descriptors_analysis_file_intermediate_bam_json
    File ss2_descriptors_analysis_file_intermediate_bai_json
    File ss2_descriptors_analysis_file_project_loom_json
    File ss2_metadata_analysis_process_file_intermediate_json
    File ss2_metadata_analysis_process_project_loom_json
    File ss2_metadata_analysis_protocol_file_intermediate_json
    File ss2_metadata_analysis_protocol_file_project_json
    File ss2_metadata_reference_file_json
    File ss2_links_json
  }

  call test_target_adapter.CreateSs2AdapterMetadata as target_adapter {
    input:
      output_bams = output_bams,
      output_bais = output_bais,
      output_loom = output_loom,
      input_ids = input_ids,
      fastq_1_uuids = fastq_1_uuids,
      fastq_2_uuids = fastq_2_uuids,
      fastq_i1_uuids = fastq_i1_uuids,
      all_libraries = all_libraries,
      all_organs = all_organs,
      all_species = all_species,
      all_project_ids = all_project_ids,
      all_project_names = all_project_names,
      workspace_bucket = workspace_bucket
  }

  # bam file descriptors come first in output_analysis_file_descriptor_objects array
  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_bam {
    input:
      test_json = target_adapter.output_analysis_file_descriptor_objects[0],
      truth_json = ss2_descriptors_analysis_file_intermediate_bam_json
  }

  # check length of array, minus project level
  # first half of array is bam descriptors, second half is bai descriptors
  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_bai {
    input:
      test_json = target_adapter.output_analysis_file_descriptor_objects[((length(target_adapter.output_analysis_file_descriptor_objects) - 1 ) / 2)],
      truth_json = ss2_descriptors_analysis_file_intermediate_bai_json
  }

  # project loom descriptor will be in last position of output_analysis_file_descriptor_objects array
  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_loom {
    input:
      test_json = target_adapter.output_analysis_file_descriptor_objects[(length(target_adapter.output_analysis_file_descriptor_objects) - 1)],
      truth_json = ss2_descriptors_analysis_file_project_loom_json
  }

  # intermediate bam analysis file is second in output_analysis_file_metadata_objects
  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_files_bam {
    input:
      test_json = target_adapter.output_analysis_file_metadata_objects[1],
      truth_json = ss2_metadata_analysis_file_intermediate_bam_json
  }

  # intermediate bai analysis files is first in output_analysis_file_metadata_objects
  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_files_bai {
    input:
      test_json = target_adapter.output_analysis_file_metadata_objects[0],
      truth_json = ss2_metadata_analysis_file_intermediate_bai_json
  }

  # project loom analysis file will be in last position of output_analysis_file_metadata_objects array
  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_files_loom {
    input:
      test_json = target_adapter.output_analysis_file_metadata_objects[(length(target_adapter.output_analysis_file_metadata_objects) - 1)],
      truth_json = ss2_metadata_analysis_file_project_loom_json
  }

  # interemediate process file is first in output_analysis_process_objects array
  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_process {
    input:
      test_json = target_adapter.output_analysis_process_objects[0],
      truth_json = ss2_metadata_analysis_process_file_intermediate_json
  }

  # project loom analysis process will be in last position of output_analysis_process_objects array
  call checker_adapter.CompareAdapterFiles as checker_adapter_ss2_project_metadata_analysis_process {
    input:
      test_json = target_adapter.output_analysis_process_objects[(length(target_adapter.output_analysis_process_objects) - 1)],
      truth_json = ss2_metadata_analysis_process_project_loom_json
  }

  # intermediate protocol file is first in output_analysis_protocol_objects array
  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_protocol {
    input:
      test_json = target_adapter.output_analysis_protocol_objects[0],
      truth_json = ss2_metadata_analysis_protocol_file_intermediate_json
  }

  # project level analysis protocol will be in last position of output_analysis_protocol_objects array
  call checker_adapter.CompareAdapterFiles as checker_adapter_ss2_project_metadata_analysis_protocol {
    input:
      test_json = target_adapter.output_analysis_protocol_objects[(length(target_adapter.output_analysis_protocol_objects) - 1)],
      truth_json = ss2_metadata_analysis_protocol_file_project_json
  }

  # only one resulting links file
  call checker_adapter.CompareAdapterFiles as checker_adapter_links {
    input:
      test_json = target_adapter.output_links_objects[0],
      truth_json = ss2_links_json
  }

  # only one resulting reference file descriptor
  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_reference {
    input:
      test_json = target_adapter.output_reference_file_descriptor_objects[0],
      truth_json = ss2_descriptors_analysis_file_intermediate_reference_json
  }

  # only one resulting reference file
  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_reference_file {
    input:
      test_json = target_adapter.output_reference_metadata_objects[0],
      truth_json = ss2_metadata_reference_file_json
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