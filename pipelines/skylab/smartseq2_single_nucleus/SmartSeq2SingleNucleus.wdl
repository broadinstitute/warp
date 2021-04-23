version 1.0

import "TrimAdapters.wdl" as TrimAdapters
import "StarAlign.wdl" as StarAlignFastq
import "Picard.wdl" as Picard
import "FeatureCounts.wdl" as CountAlignments


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

    String stranded
    String input_id
    String output_name

    File adapter_list
    File fastq1
    File? fastq2
    Boolean paired_end
  }
  
  # version of this pipeline
  String pipeline_version = "1.0.0"

  parameter_meta {
    stranded: "Library strand information example values: FR RF NONE"
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

  output {
    # version of this pipeline
    String pipeline_version_out = pipeline_version

    # quality control outputs
    # data outputs
    # loom
  }
}
