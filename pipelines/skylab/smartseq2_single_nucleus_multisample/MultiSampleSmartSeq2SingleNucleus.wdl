version 1.0

import "../../../pipelines/skylab/smartseq2_single_nucleus/SmartSeq2SingleNucleus.wdl" as single_nucleus_run
import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils
import "../../../tasks/skylab/CheckInputs.wdl" as CheckInputs
       
workflow MultiSampleSmartSeq2SingleNucleus {
  meta {
    description: "The MultiSampleSmartSeq2SingleNucleus pipeline runs multiple snSS2 samples in a single pipeline invocation"
    allowNestedInputs: true
  }

  input {
      # reference genome fasta
      File genome_ref_fasta

      # Reference index information
      File tar_star_reference
      # annotation file
      File annotations_gtf
      # adapter list file
      File adapter_list

      # Sample information
      Array[String] input_ids
      Array[String]? input_names
      Array[String] fastq1_input_files
      Array[String] fastq2_input_files 
      String batch_id
      String? batch_name
      Array[String]? project_id
      Array[String]? project_name
      Array[String]? library
      Array[String]? species
      Array[String]? organ
      String? input_name_metadata_field
      String? input_id_metadata_field
  }
  # Version of this pipeline
  String pipeline_version = "1.0.0"

  if (false) {
     String? none = "None"
  }

  # Parameter metadata information
  parameter_meta {
    genome_ref_fasta: "Genome reference in fasta format"
    tar_star_reference: "STAR reference index tar file"
    annotations_gtf: "gtf containing annotations for gene tagging (must match star reference)"
    input_ids: "Array of input ids"
    input_names: "Array of input names"
    input_id_metadata_field: "String that describes the metadata field containing the input_ids"
    input_name_metadata_field: "String that describes the metadata field containing the input_names"
    fastq1_input_files: "Array of fastq1 files; order must match the order in input_id."
    fastq2_input_files: "Array of fastq2 files for paired end runs; order must match fastq1_input_files and input_id."
    batch_id: " Identifier for the batch"
  }

  # Check that all input arrays are the same length
  call CheckInputs.checkInputArrays as checkArrays{
      input:
         input_ids = input_ids,
         input_names = input_names,
         fastq1_input_files = fastq1_input_files,
         fastq2_input_files = fastq2_input_files
  }

  ### Execution starts here ###
  scatter(idx in range(length(input_ids))) {
      call single_nucleus_run.SmartSeq2SingleNucleus as sn_pe {
        input:
           genome_ref_fasta = genome_ref_fasta,
           star_reference = tar_star_reference,
           annotations_gtf =  annotations_gtf,
           adapter_list = adapter_list,
           fastq1 = fastq1_input_files[idx],
           fastq2 = fastq2_input_files[idx],
           input_id = input_ids[idx],
           output_name = input_ids[idx],
           input_name_metadata_field = input_name_metadata_field,
           input_id_metadata_field = input_id_metadata_field,
           input_name = if defined(input_names) then select_first([input_names])[idx] else none
      }
  }


  ### Aggregate the Loom Files Directly ###
  call LoomUtils.AggregateSmartSeq2Loom as AggregateLoom {
    input:
      loom_input = sn_pe.loom_output_file,
      batch_id = batch_id,
      batch_name = batch_name,
      project_id = if defined(project_id) then select_first([project_id])[0] else none,
      project_name = if defined(project_name) then select_first([project_name])[0] else none,
      library = if defined(library) then select_first([library])[0] else none,
      species = if defined(species) then select_first([species])[0] else none,
      organ = if defined(organ) then select_first([organ])[0] else none,
      pipeline_version = "MultiSampleSmartSeq2SingleNucleus_v~{pipeline_version}"
  }


  ### Pipeline output ###
  output {
    # loom output, exon/intron count tsv files and the aligned bam files
    File loom_output = AggregateLoom.loom_output_file
    Array[File] exon_intron_count_files = sn_pe.exon_intron_counts 
    Array[File] bam_files = sn_pe.aligned_bam
  }
}
