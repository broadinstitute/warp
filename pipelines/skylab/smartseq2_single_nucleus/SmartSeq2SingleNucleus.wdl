version 1.0

import "https://raw.githubusercontent.com/broadinstitute/warp/snSS2_first_wdls/tasks/skylab/TrimAdapters.wdl" as TrimAdapters
import "https://raw.githubusercontent.com/broadinstitute/warp/snSS2_first_wdls/tasks/skylab/StarAlign.wdl" as StarAlignFastq
import "https://raw.githubusercontent.com/broadinstitute/warp/snSS2_first_wdls/tasks/skylab/Picard.wdl" as Picard
import "https://raw.githubusercontent.com/broadinstitute/warp/snSS2_first_wdls/tasks/skylab/FeatureCounts.wdl" as CountAlignments
import "https://raw.githubusercontent.com/broadinstitute/warp/snSS2_first_wdls/tasks/skylab/LoomUtils.wdl" as LoomUtils


workflow SmartSeq2SingleNucleus {
  meta {
    description: "Process SmartSeq2 snRNA-Seq data, include reads alignment, QC metrics collection, and gene expression quantitication"
  }

  input {
    # load annotation
    File genome_ref_fasta
    # load index
    # ref index name

    File tar_star_reference
    File annotations_gtf

    # samples
    String input_id

    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field

    String output_name

    File adapter_list
    File fastq1
    File? fastq2
    Boolean paired_end
  }
  # version of this pipeline
  String pipeline_version = "1.0.0"

  parameter_meta {
    input_id: "Sample name or cell_names"
    output_name: "Output name, can include path"
    fastq1: "R1 in paired end reads"
    fastq2: "R2 in paired end reads"
    paired_end: "Boolean flag denoting if the sample is paired end or not"
    tar_star_reference: "star genome reference"
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
        tar_star_reference = tar_star_reference
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

  Array[File] picard_row_outputs = [
     CollectMultipleMetrics.alignment_summary_metrics,
     RemoveDuplicatesFromBam.dedup_metrics,
     CollectMultipleMetrics.gc_bias_summary_metrics
  ]

  call LoomUtils.SingleNucleiSmartSeq2LoomOutput as SingleNucleiSmartSeq2LoomOutput {
       input:
         input_id=input_id,
         input_name = input_name,
         pipeline_version = "SmartSeq2SingleNucleus_v~{pipeline_version}",
         input_id_metadata_field = input_id_metadata_field,
         input_name_metadata_field = input_name_metadata_field,
         smartseq_qc_files = picard_row_outputs,
         introns_counts = CountAlignments.intron_counts_out,
         exons_counts = CountAlignments.exon_counts_out,
         annotation_introns_added_gtf = annotations_gtf
  }

  output {
    # version of this pipeline
    String pipeline_version_out = pipeline_version
    # duplicate removed BAM
    File aligned_bam = RemoveDuplicatesFromBam.output_bam
    # data outputs
    File exon_intron_counts=SingleNucleiSmartSeq2LoomOutput.exon_intron_counts
    # loom
    File loom_output_files = SingleNucleiSmartSeq2LoomOutput.loom_output
  }
}