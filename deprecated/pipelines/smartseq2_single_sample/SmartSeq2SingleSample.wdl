version 1.0

# SmartSeq2SingleSample is now deprecated 2025-03-06

import "../../../tasks/wdl/HISAT2.wdl" as HISAT2
import "../../../tasks/wdl/Picard.wdl" as Picard
import "../../../tasks/wdl/RSEM.wdl" as RSEM
import "../../../tasks/wdl/GroupMetricsOutputs.wdl" as GroupQCs
import "../../../tasks/wdl/LoomUtils.wdl" as LoomUtils

workflow SmartSeq2SingleSample {
  meta {
    description: "Process SmartSeq2 scRNA-Seq data, include reads alignment, QC metrics collection, and gene expression quantitication"
  }

  input {
    # load annotation
    File genome_ref_fasta
    File rrna_intervals
    File gene_ref_flat
    # load index
    File hisat2_ref_index
    File hisat2_ref_trans_index
    File rsem_ref_index
    # ref index name
    String hisat2_ref_name
    String hisat2_ref_trans_name
    # samples
    String stranded
    String input_id
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field
    String output_name
    File fastq1
    File? fastq2
    Boolean paired_end
  }
  
  # version of this pipeline
  String pipeline_version = "5.1.22"

  parameter_meta {
    genome_ref_fasta: "Genome reference in fasta format"
    rrna_intervals: "rRNA interval file required by Picard"
    gene_ref_flat: "Gene refflat file required by Picard"
    hisat2_ref_index: "HISAT2 reference index file in tarball"
    hisat2_ref_trans_index: "HISAT2 transcriptome index file in tarball"
    rsem_ref_index: "RSEM reference index file in tarball"
    hisat2_ref_name: "HISAT2 reference index name"
    hisat2_ref_trans_name: "HISAT2 transcriptome index file name"
    stranded: "Library strand information example values: FR RF NONE"
    input_id: "Sample name or cell_names"
    input_id_metadata_field: "String that describes the metadata field containing the input_id"
    input_name: "User provided sample name or cell_names"
    input_name_metadata_field: "String that describes the metadata field containing input_name"
    output_name: "Output name, can include path"
    fastq1: "R1 in paired end reads"
    fastq2: "R2 in paired end reads"
    paired_end: "Boolean flag denoting if the sample is paired end or not"
  }

  String quality_control_output_basename = output_name + "_qc"

   if( paired_end ) {
     call HISAT2.HISAT2PairedEnd {
       input:
         hisat2_ref = hisat2_ref_index,
         fastq1 = fastq1,
         fastq2 = select_first([fastq2]),
         ref_name = hisat2_ref_name,
         input_id = input_id,
         output_basename = quality_control_output_basename,
     }
  }
  if( !paired_end ) {
     call HISAT2.HISAT2SingleEnd {
       input:
         hisat2_ref = hisat2_ref_index,
         fastq = fastq1,
         ref_name = hisat2_ref_name,
         input_id = input_id,
         output_basename = quality_control_output_basename,
     }
  }

  File HISAT2_output_bam = select_first([ HISAT2PairedEnd.output_bam, HISAT2SingleEnd.output_bam] )
  File HISAT2_bam_index = select_first([ HISAT2PairedEnd.bam_index, HISAT2SingleEnd.bam_index] )
  File HISAT2_log_file = select_first([ HISAT2PairedEnd.log_file, HISAT2SingleEnd.log_file] )

  call Picard.CollectMultipleMetrics {
    input:
      aligned_bam = HISAT2_output_bam,
      genome_ref_fasta = genome_ref_fasta,
      output_basename = quality_control_output_basename,
  }

  call Picard.CollectRnaMetrics {
    input:
      aligned_bam = HISAT2_output_bam,
      ref_flat = gene_ref_flat,
      rrna_intervals = rrna_intervals,
      output_basename = quality_control_output_basename,
      stranded = stranded,
  }

  call Picard.CollectDuplicationMetrics {
    input:
      aligned_bam = HISAT2_output_bam,
      output_basename = quality_control_output_basename,
  }

  String data_output_basename = output_name + "_rsem"

  if( paired_end ) {
      call HISAT2.HISAT2RSEM as HISAT2Transcriptome {
        input:
          hisat2_ref = hisat2_ref_trans_index,
          fastq1 = fastq1,
          fastq2 = fastq2,
          ref_name = hisat2_ref_trans_name,
          input_id = input_id,
          output_basename = data_output_basename,
      }
  }

  if( !paired_end ) {
      call HISAT2.HISAT2RSEMSingleEnd as HISAT2SingleEndTranscriptome {
        input:
          hisat2_ref = hisat2_ref_trans_index,
          fastq = fastq1,
          ref_name = hisat2_ref_trans_name,
          input_id = input_id,
          output_basename = data_output_basename,
      }
  }

  File HISAT2RSEM_output_bam = select_first([ HISAT2Transcriptome.output_bam, HISAT2SingleEndTranscriptome.output_bam] )
  File HISAT2RSEM_log_file = select_first([ HISAT2Transcriptome.log_file, HISAT2SingleEndTranscriptome.log_file] )

  call RSEM.RSEMExpression {
    input:
      trans_aligned_bam = HISAT2RSEM_output_bam,
      rsem_genome = rsem_ref_index,
      output_basename = data_output_basename,
      is_paired = paired_end
  }

  Array[File] picard_row_outputs = [CollectMultipleMetrics.alignment_summary_metrics,CollectDuplicationMetrics.dedup_metrics,CollectRnaMetrics.rna_metrics,CollectMultipleMetrics.gc_bias_summary_metrics]

  # This output only exists for PE and select_first fails if array is empty
  if ( length(CollectMultipleMetrics.insert_size_metrics) > 0 ) {
    File? picard_row_optional_outputs = select_first(CollectMultipleMetrics.insert_size_metrics)
  }

  Array[File] picard_table_outputs = [
    CollectMultipleMetrics.base_call_dist_metrics,
    CollectMultipleMetrics.gc_bias_detail_metrics,
    CollectMultipleMetrics.pre_adapter_details_metrics,
    CollectMultipleMetrics.pre_adapter_summary_metrics,
    CollectMultipleMetrics.bait_bias_detail_metrics,
    CollectMultipleMetrics.error_summary_metrics,
  ]

  call GroupQCs.GroupQCOutputs {
   input:
      picard_row_outputs = picard_row_outputs,
      picard_row_optional_outputs = select_all(CollectMultipleMetrics.insert_size_metrics),
      picard_table_outputs = picard_table_outputs,
      hisat2_stats = HISAT2_log_file,
      hisat2_trans_stats = HISAT2RSEM_log_file,
      rsem_stats = RSEMExpression.rsem_cnt,
      output_name = output_name
  }

  call LoomUtils.SmartSeq2LoomOutput {
    input:
      rsem_gene_results = RSEMExpression.rsem_gene,
      smartseq_qc_files = GroupQCOutputs.group_files,
      input_id=input_id,
      input_name = input_name,
      pipeline_version = "SmartSeq2SingleSample_v~{pipeline_version}",
      input_id_metadata_field = input_id_metadata_field,
      input_name_metadata_field = input_name_metadata_field
  }

  output {
    # version of this pipeline
    String pipeline_version_out = pipeline_version

    # quality control outputs
    File aligned_bam = HISAT2_output_bam
    File bam_index = HISAT2_bam_index
    File? insert_size_metrics =  picard_row_optional_outputs
    File quality_distribution_metrics = CollectMultipleMetrics.quality_distribution_metrics
    File quality_by_cycle_metrics = CollectMultipleMetrics.quality_by_cycle_metrics
    File bait_bias_summary_metrics = CollectMultipleMetrics.bait_bias_summary_metrics
    File rna_metrics = CollectRnaMetrics.rna_metrics # check this
    Array[File] group_results = GroupQCOutputs.group_files # check this

    # data outputs
    File aligned_transcriptome_bam = HISAT2RSEM_output_bam
    File rsem_gene_results = RSEMExpression.rsem_gene
    File rsem_isoform_results = RSEMExpression.rsem_isoform

    # loom
    File loom_output_files = SmartSeq2LoomOutput.loom_output
  }
}
