version 1.0

import "../../../tasks/skylab/TrimAdapters.wdl" as TrimAdapters
import "../../../tasks/skylab/StarAlign.wdl" as StarAlignFastq
import "../../../tasks/skylab/FeatureCounts.wdl" as CountAlignments
import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils
import "../../../tasks/skylab/Picard.wdl" as Picard


workflow SmartSeq2SingleNucleus {
  meta {
    description: "Process SmartSeq2 snRNA-Seq data, including read alignments, QC metrics collection, and gene expression quantification, for intronic and exonic regions"
  }

  input {
    # reference genome FASTA file 
    File genome_ref_fasta

    # STAR ref index name
    File star_reference
    # annotation file 
    File annotations_gtf

    # sample id
    String input_id

    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field

    String output_name

    File adapter_list
    # at this point only paired_end reads are supported
    File fastq1
    File fastq2
  }
  # version of this pipeline
  String pipeline_version = "1.0.0"

  parameter_meta {
    input_id: "Sample name or cell_names"
    output_name: "Output name, can include path"
    fastq1: "R1 in paired end reads"
    fastq2: "R2 in paired end reads"
    star_reference: "star genome reference in the form of a tar file"
    annotations_gtf: "gtf containing annotations for gene tagging (must match star reference)"
  }

  String quality_control_output_basename = output_name + "_qc"

  call TrimAdapters.TrimAdapters as TrimAdapters {
       input:
         fastq1 = fastq1,
         fastq2 = select_first([fastq2]),
         adapter_list = adapter_list
   }

   call StarAlignFastq.StarAlignFastqPairedEnd as StarAlign {
      input:
        fastq1 = TrimAdapters.trimmed_fastq1,
        fastq2 = TrimAdapters.trimmed_fastq2,
        tar_star_reference = star_reference
   }

  call Picard.RemoveDuplicatesFromBam as RemoveDuplicatesFromBam {
    input:
      input_id = input_id,
      aligned_bam = StarAlign.output_bam,
      output_basename = quality_control_output_basename,
  }

  call Picard.CollectMultipleMetrics {
    input:
      aligned_bam = RemoveDuplicatesFromBam.output_bam,
      genome_ref_fasta = genome_ref_fasta,
      output_basename = quality_control_output_basename,
  }

  call CountAlignments.CountAlignments as CountAlignments {
       input:
         input_bam = RemoveDuplicatesFromBam.output_bam,
         annotation_gtf = annotations_gtf
   }

  Array[File] smartseq_qc_files = [
     CollectMultipleMetrics.alignment_summary_metrics,
     RemoveDuplicatesFromBam.dedup_metrics,
     CollectMultipleMetrics.gc_bias_summary_metrics
  ]

  call LoomUtils.SingleNucleiSmartSeq2LoomOutput as SingleNucleiSmartSeq2LoomOutput {
       input:
         input_id = input_id,
         input_name = input_name,
         pipeline_version = "SmartSeq2SingleNucleus_v~{pipeline_version}",
         input_id_metadata_field = input_id_metadata_field,
         input_name_metadata_field = input_name_metadata_field,
         smartseq_qc_files = smartseq_qc_files,
         introns_counts = CountAlignments.intron_counts_out,
         exons_counts = CountAlignments.exon_counts_out,
         annotation_introns_added_gtf = annotations_gtf
  }

  output {
    # version of this pipeline
    String pipeline_version_out = pipeline_version
    # deduplicated BAM
    File aligned_bam = RemoveDuplicatesFromBam.output_bam
    # data outputs
    File exon_intron_counts=SingleNucleiSmartSeq2LoomOutput.exon_intron_counts
    # this loom contains the exonic and intronic counts for a sample corresponding to individual genes, and metadata
    File loom_output_file = SingleNucleiSmartSeq2LoomOutput.loom_output
  }
}
