version 1.0
# MultiSampleSmartSeq2 is now deprecated 2025-03-06

import "../../../pipelines/wdl/smartseq2_single_sample/SmartSeq2SingleSample.wdl" as single_cell_run
import "../../../tasks/wdl/LoomUtils.wdl" as LoomUtils
       
workflow MultiSampleSmartSeq2 {
  meta {
    description: "The MultiSampleSmartSeq2 pipeline runs multiple SS2 samples in a single pipeline invocation"
    allowNestedInputs: true
  }

  input {
      # Gene Annotation
      File genome_ref_fasta
      File rrna_intervals
      File gene_ref_flat

      # Reference index information
      File hisat2_ref_name
      File hisat2_ref_trans_name
      File hisat2_ref_index
      File hisat2_ref_trans_index
      File rsem_ref_index

      # Sample information
      String stranded
      Array[String] input_ids
      Array[String]? input_names
      Array[String] fastq1_input_files
      Array[String] fastq2_input_files = []
      String batch_id
      String? batch_name
      Array[String]? project_id
      Array[String]? project_name
      Array[String]? library
      Array[String]? species
      Array[String]? organ
      String? input_name_metadata_field
      String? input_id_metadata_field
      Boolean paired_end
  }
  # Version of this pipeline
  String pipeline_version = "2.2.23"

  if (false) {
     String? none = "None"
  }

  # Parameter metadata information
  parameter_meta {
    genome_ref_fasta: "Genome reference in fasta format"
    rrna_intervals: "rRNA interval file required by Picard"
    gene_ref_flat: "Gene refflat file required by Picard"
    hisat2_ref_name: "HISAT2 reference index name"
    hisat2_ref_trans_name: "HISAT2 transcriptome index file name"
    hisat2_ref_index: "HISAT2 reference index file in tarball"
    hisat2_ref_trans_index: "HISAT2 transcriptome index file in tarball"
    rsem_ref_index: "RSEM reference index file in tarball"
    stranded: "Library strand information example values: FR RF NONE"
    input_ids: "Array of input ids"
    input_names: "Array of input names"
    input_id_metadata_field: "String that describes the metadata field containing the input_ids"
    input_name_metadata_field: "String that describes the metadata field containing the input_names"
    fastq1_input_files: "Array of fastq1 files; order must match the order in input_id."
    fastq2_input_files: "Array of fastq2 files for paired end runs; order must match fastq1_input_files and input_id."
    batch_id: " Identifier for the batch"
    paired_end: "Is the sample paired end or not"
  }

  # Check that all input arrays are the same length
  call checkInputArrays as checkArrays{
      input:
         paired_end = paired_end,
         input_ids = input_ids,
         input_names = input_names,
         fastq1_input_files = fastq1_input_files,
         fastq2_input_files = fastq2_input_files
  }

  ### Execution starts here ###
  if (paired_end) {
    scatter(idx in range(length(input_ids))) {
      call single_cell_run.SmartSeq2SingleSample as sc_pe {
        input:
          fastq1 = fastq1_input_files[idx],
          fastq2 = fastq2_input_files[idx],
          stranded = stranded,
          genome_ref_fasta = genome_ref_fasta,
          rrna_intervals = rrna_intervals,
          gene_ref_flat = gene_ref_flat,
          hisat2_ref_index = hisat2_ref_index,
          hisat2_ref_name = hisat2_ref_name,
          hisat2_ref_trans_index = hisat2_ref_trans_index,
          hisat2_ref_trans_name = hisat2_ref_trans_name,
          rsem_ref_index = rsem_ref_index,
          input_id = input_ids[idx],
          output_name = input_ids[idx],
          paired_end = paired_end,
          input_name_metadata_field = input_name_metadata_field,
          input_id_metadata_field = input_id_metadata_field,
          input_name = if defined(input_names) then select_first([input_names])[idx] else none
      }
    }
  }
  if (!paired_end) {
    scatter(idx in range(length(input_ids))) {
      call single_cell_run.SmartSeq2SingleSample as sc_se {
        input:
          fastq1 = fastq1_input_files[idx],
          stranded = stranded,
          genome_ref_fasta = genome_ref_fasta,
          rrna_intervals = rrna_intervals,
          gene_ref_flat = gene_ref_flat,
          hisat2_ref_index = hisat2_ref_index,
          hisat2_ref_name = hisat2_ref_name,
          hisat2_ref_trans_index = hisat2_ref_trans_index,
          hisat2_ref_trans_name = hisat2_ref_trans_name,
          rsem_ref_index = rsem_ref_index,
          input_id = input_ids[idx],
          output_name = input_ids[idx],
          paired_end = paired_end,
          input_name_metadata_field = input_name_metadata_field,
          input_id_metadata_field = input_id_metadata_field,
          input_name = if defined(input_names) then select_first([input_names])[idx] else none

      }
    }
  }

  Array[File] loom_output_files = select_first([sc_pe.loom_output_files, sc_se.loom_output_files])
  Array[File] bam_files_intermediate = select_first([sc_pe.aligned_bam, sc_se.aligned_bam])
  Array[File] bam_index_files_intermediate = select_first([sc_pe.bam_index, sc_se.bam_index])

  ### Aggregate the Loom Files Directly ###
  call LoomUtils.AggregateSmartSeq2Loom as AggregateLoom {
    input:
      loom_input = loom_output_files,
      batch_id = batch_id,
      batch_name = batch_name,
      project_id = if defined(project_id) then select_first([project_id])[0] else none,
      project_name = if defined(project_name) then select_first([project_name])[0] else none,
      library = if defined(library) then select_first([library])[0] else none,
      species = if defined(species) then select_first([species])[0] else none,
      organ = if defined(organ) then select_first([organ])[0] else none,
      pipeline_version = "MultiSampleSmartSeq2_v~{pipeline_version}"
  }


  ### Pipeline output ###
  output {
    # Bam files and their indexes
    Array[File] bam_files = bam_files_intermediate
    Array[File] bam_index_files = bam_index_files_intermediate
    File loom_output = AggregateLoom.loom_output_file
  }
}

task checkInputArrays {
  input {
    Boolean paired_end
    Array[String] input_ids
    Array[String]? input_names
    Array[String] fastq1_input_files
    Array[String] fastq2_input_files
  }
  Int len_input_ids = length(input_ids)
  Int len_fastq1_input_files = length(fastq1_input_files)
  Int len_fastq2_input_files = length(fastq2_input_files)
  Int len_input_names = if defined(input_names) then length(select_first([input_names])) else 0

  meta {
    description: "checks input arrays to ensure that all arrays are the same length"
  }

  command {
    set -e

    if [[ ~{len_input_ids} !=  ~{len_fastq1_input_files} ]]
      then
      echo "ERROR: Different number of arguments for input_id and fastq1 files"
      exit 1;
    fi

    if [[ ~{len_input_names} != 0  && ~{len_input_ids} !=  ~{len_input_names} ]]
        then
        echo "ERROR: Different number of arguments for input_name and input_id"
        exit 1;
    fi

    if  ~{paired_end} && [[ ~{len_fastq2_input_files} != ~{len_input_ids} ]]
      then
      echo "ERROR: Different number of arguments for sample names and fastq1 files"
      exit 1;
    fi
    exit 0;
  }

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk 1 HDD"
  }

}
