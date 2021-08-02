version 1.0

import "../../../../projects/hca_mvp/getTerraMetadata.wdl" as target_adapter
import "../../../../tests/skylab/hca_adapter/pr/ValidateHcaAdapter.wdl" as checker_adapter
import "../../../../projects/hca_mvp/OptimusPostProcessing.wdl" as target_OptimusPostProcessing

# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestHcaAdapter {
  input {

    # hca adapter inputs
    Boolean add_md5s
    String analysis_file_version
    String analysis_process_schema_version
    String analysis_protocol_schema_version
    String assembly_type
    String cromwell_url
    Array[String] fastq_1_uuids
    Array[String] fastq_2_uuids
    String file_descriptor_schema_version
    File genome_ref_fasta
    String genus_species
    Array[String] input_uuids
    Array[String] library
    String links_schema_version
    String method
    String ncbi_taxon_id
    Array[String] organ
    Array[File] outputs
    String pipeline_tools_version
    String pipeline_version
    Array[String] projects
    String reference_file_schema_version
    String reference_type
    String reference_version
    String run_type
    String schema_url
    Array[String] species
    String staging_area

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

    #optimus post processing inputs
    Array[File] library_looms
    Array[File] analysis_file_jsons
    Array[File] links_jsons
    Array[String] post_processing_library
    Array[String] species
    Array[String] organ
    String project_id
    String project_name
    String output_basename
    String post_processing_staging_area = "gs://fc-c307d7b3-8386-40a1-b32c-73b9e16e0103/"
    String version_timestamp = "2021-05-24T12:00:00.000000Z"


  }


  call target_adapter.submit as target_adapter {
    input:
    add_md5s=add_md5s,
    analysis_file_version=analysis_file_version,
    analysis_process_schema_version=analysis_process_schema_version,
    analysis_protocol_schema_version=analysis_protocol_schema_version,
    assembly_type=assembly_type,
    cromwell_url=cromwell_url,
    fastq_1_uuids=fastq_1_uuids,
    fastq_2_uuids=fastq_2_uuids,
    file_descriptor_schema_version=file_descriptor_schema_version,
    genome_ref_fasta=genome_ref_fasta,
    genus_species=genus_species,
    input_uuids=input_uuids,
    library=library,
    links_schema_version=links_schema_version,
    method=method,
    ncbi_taxon_id=ncbi_taxon_id,
    organ=organ,
    outputs=outputs,
    pipeline_tools_version=pipeline_tools_version,
    pipeline_version=pipeline_version,
    projects=projects,
    reference_file_schema_version=reference_file_schema_version,
    reference_type=reference_type,
    reference_version=reference_version,
    run_type=run_type,
    schema_url=schema_url,
    species=species,
    staging_area=staging_area
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_loom {
   input:
     test_json=target_adapter.analysis_file_descriptor[0],
     truth_json=optimus_descriptors_analysis_file_intermediate_loom_json
   }

  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_bam {
   input:
     test_json=target_adapter.analysis_file_descriptor[1],
     truth_json=optimus_descriptors_analysis_file_intermediate_bam_json
   }

  call checker_adapter.CompareAdapterFiles as checker_adapter_descriptors_reference{
    input:
      test_json=target_adapter.reference_genome_descriptor,
      truth_json=optimus_descriptors_analysis_file_intermediate_reference_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_links {
    input:
     test_json=target_adapter.links,
     truth_json=optimus_links_intermediate_loom_json,
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_files_loom {
    input:
     test_json=target_adapter.analysis_file[0],
     truth_json=optimus_metadata_analysis_file_intermediate_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_files_bam {
    input:
      test_json=target_adapter.analysis_file[1],
      truth_json=optimus_metadata_analysis_file_intermediate_bam_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_process {
    input:
      test_json=target_adapter.analysis_process,
      truth_json=optimus_metadata_analysis_process_file_intermediate_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_reference_file {
      input:
        test_json=target_adapter.reference_genome_reference_file,
        truth_json=optimus_metadata_reference_file_intermediate_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_metadata_analysis_protocol {
    input:
      test_json=target_adapter.analysis_protocol,
      truth_json=optimus_metadata_analysis_protocol_file_intermediate_json

  }

  call target_OptimusPostProcessing.OptimusPostProcessing as target_OptimusPostProcessing {
    input:
      library_looms=library_looms,
      analysis_file_jsons=analysis_file_jsons,
      links_jsons=links_jsons,
      post_processing_library=post_processing_library,
      species=species,
      organ=organ,
      project_id=project_id,
      project_name=project_name,
      output_basename=output_basename,
      post_processing_staging_area=post_processing_staging_area,
	  version_timestamp=version_timestamp
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_descriptors {
       input:
         test_json=target_OptimusPostProcessing.json_adapter_files[3],
         truth_json=optimus_descriptors_analysis_file_project_loom_json
  }

  call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_links {
       input:
         test_json=target_OptimusPostProcessing.json_adapter_files[4],
         truth_json=optimus_links_project_loom_json
  }

 call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_metadata {
      input:
        test_json=target_OptimusPostProcessing.json_adapter_files[0],
        truth_json=optimus_metadata_analysis_file_project_loom_json
 }

 call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_metadata_analysis_process {
      input:
        test_json=target_OptimusPostProcessing.json_adapter_files[1],
        truth_json=optimus_metadata_analysis_process_project_loom_json
 }

 call checker_adapter.CompareAdapterFiles as checker_adapter_optimus_project_metadata_analysis_protocol {
       input:
         test_json=target_OptimusPostProcessing.json_adapter_files[2],
         truth_json=optimus_metadata_analysis_protocol_project_loom_json
 }

   output {
    Array[File] analysis_file = target_adapter.analysis_file
    File analysis_process = target_adapter.analysis_process
    File analysis_protocol = target_adapter.analysis_protocol
    Array[File] analysis_output = target_adapter.analysis_output
    File reference_genome = target_adapter.reference_genome
    File reference_genome_reference_file = target_adapter.reference_genome_reference_file
    File reference_genome_descriptor = target_adapter.reference_genome_descriptor
    Array[File] analysis_file_descriptor = target_adapter.analysis_file_descriptor
    File links = target_adapter.links
   }
}