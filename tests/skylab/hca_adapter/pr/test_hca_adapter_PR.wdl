version 1.0

import "https://raw.githubusercontent.com/broadinstitute/warp/np-hca-adapter-test/pipelines/skylab/optimus/Optimus.wdl" as target_optimus
import "https://raw.githubusercontent.com/broadinstitute/warp/np-hca-adapter-test/beta-pipelines/skylab/hca_adapter/getTerraMetadata.wdl" as target_adapter
import "https://raw.githubusercontent.com/broadinstitute/warp/np-hca-adapter-test/tests/skylab/hca_adapter/pr/ValidateHcaAdapter.wdl" as checker_adapter

# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestHcaAdapter {
  input {

    # Optimus inputs
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq

    File whitelist  # 10x genomics cell barcode whitelist for 10x V2
    File tar_star_reference  # star reference
    File annotations_gtf  # gtf containing annotations for gene tagging
    File ref_genome_fasta  # genome fasta file
    String input_id  # name of sample matching this file, inserted into read group header
    String chemistry

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
  }

  call target_optimus.Optimus as target_optimus {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      i1_fastq = i1_fastq,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      input_id = input_id,
      emptydrops_lower = 1,
      chemistry = chemistry
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
    outputs=[target_optimus.loom_output_file, target_optimus.bam],
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

  call checker_adapter.ValidateOptimusDescriptorAnalysisFiles as checker_adapter_descriptors {
     input:
       optimus_descriptors_analysis_file_intermediate_loom_json=target_adapter.analysis_file_descriptor[1],
       expected_optimus_descriptors_analysis_file_intermediate_loom_json_hash='007e56d8bf4f785ee37e9be2325a534c',
       optimus_descriptors_analysis_file_intermediate_bam_json=target_adapter.analysis_file_descriptor[0],
       expected_optimus_descriptors_analysis_file_intermediate_bam_json_hash='8132649e8263f7d59d466291761d9ff0',
       optimus_descriptors_analysis_file_intermediate_reference_json=target_adapter.reference_genome_descriptor,
       expected_optimus_descriptors_analysis_file_intermediate_reference_json_hash='37041798cf2bb4ea4242572cbc7cda96'
   }

  call checker_adapter.ValidateOptimusLinksFiles as checker_adapter_links {
    input:
     optimus_links_intermediate_loom_json=target_adapter.links,
     expected_optimus_links_intermediate_loom_json='eebddd1aa0877ff67d446de352409dc6'
  }

  call checker_adapter.ValidateOptimusMetadataAnalysisFiles as checker_adapter_metadata_analysis_files {
    input:
     optimus_metadata_analysis_file_intermediate_bam_json=target_adapter.analysis_file[0],
     expected_optimus_metadata_analysis_file_intermediate_bam_json_hash='a259e8a60e823689550e9dc8e6e0ce04',
     optimus_metadata_analysis_file_intermediate_loom_json=target_adapter.analysis_file[1],
     expected_optimus_metadata_analysis_file_intermediate_loom_json_hash='4361d3bd65023b0dd42f35b549c0b087'
  }

  call checker_adapter.ValidateOptimusMetadataAnalysisProcessFiles as checker_adapter_metadata_analysis_process {
    input:
      optimus_metadata_analysis_process_file_intermediate_json=target_adapter.analysis_process,
      expected_optimus_metadata_analysis_process_file_intermediate_json_hash='650ae50f5f1cf886f275ada3c69ce4c2'
  }

  call checker_adapter.ValidateOptimusMetadataAnalysisProtocolFiles as checker_adapter_metadata_analysis_protocol {
    input:
      optimus_metadata_analysis_protocol_file_intermediate_json=target_adapter.analysis_protocol,
      expected_optimus_metadata_analysis_protocol_file_intermediate_json_hash='ea95eb817c461594412b56db970917a2'
  }

  call checker_adapter.ValidateOptimusMetadataReferenceFiles as checker_adapter_metadata_reference_file {
      input:
        optimus_metadata_reference_file_intermediate_json=target_adapter.reference_genome_reference_file,
        expected_optimus_metadata_reference_file_intermediate_json_hash='fakehash'
    }

}