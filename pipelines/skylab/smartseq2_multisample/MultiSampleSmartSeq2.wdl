version 1.0

import "../../../pipelines/skylab/smartseq2_single_sample/SmartSeq2SingleSample.wdl" as single_cell_run
import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils
       
workflow MultiSampleSmartSeq2 {
  meta {
    description: "The MultiSampleSmartSeq2 pipeline runs multiple SS2 samples in a single pipeline invocation"
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
      Array[LoomMetadata] input_fields
      Array[String] fastq1_input_files
      Array[String] fastq2_input_files = []
      String batch_id
      String? batch_name
      Boolean paired_end
  }
  # Version of this pipeline
  String pipeline_version = "2.1.0"

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
    fastq1_input_files: "Array of fastq1 files; order must match the order in input_id."
    fastq2_input_files: "Array of fastq2 files for paired end runs; order must match fastq1_input_files and input_id."
    batch_id: " Identifier for the batch"
    paired_end: "Is the sample paired end or not"
  }

  # Check that all input arrays are the same length
  call checkInputArrays as checkArrays{
      input:
         paired_end = paired_end,
         input_fields = input_fields,
         fastq1_input_files = fastq1_input_files,
         fastq2_input_files = fastq2_input_files
  }

  ### Execution starts here ###
  if (paired_end) {
    scatter(idx in range(length(input_fields))) {
      call single_cell_run.SmartSeq2SingleCell as sc_pe {
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
          input_fields = input_fields[idx],
          paired_end = paired_end
      }
    }
  }
  if (!paired_end) {
    scatter(idx in range(length(input_fields))) {
      call single_cell_run.SmartSeq2SingleCell as sc_se {
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
          input_fields = input_fields[idx],
          paired_end = paired_end

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
    Array[LoomMetadata] input_fields
    Array[String] fastq1_input_files
    Array[String] fastq2_input_files
  }
  Int len_input_ids = length(input_fields)
  Int len_fastq1_input_files = length(fastq1_input_files)
  Int len_fastq2_input_files = length(fastq2_input_files)

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
