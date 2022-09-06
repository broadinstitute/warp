version 1.0

import "../../../tasks/skylab/CheckInputs.wdl" as CheckInputs
import "../../../tasks/skylab/TrimAdapters.wdl" as TrimAdapters
import "../../../tasks/skylab/StarAlign.wdl" as StarAlign
import "../../../tasks/skylab/Picard.wdl" as Picard
import "../../../tasks/skylab/FeatureCounts.wdl" as CountAlignments
import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils

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
  String pipeline_version = "1.2.13"

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
         fastq2_input_files = fastq2_input_files,
         paired_end = true
  }
  call TrimAdapters.TrimAdapters as TrimAdapters {
       input:
         input_ids = input_ids,
         fastq1_input_files = fastq1_input_files,
         fastq2_input_files = fastq2_input_files,
         adapter_list = adapter_list
   }

   call StarAlign.StarAlignFastqMultisample as StarAlign {
      input:
        input_ids = input_ids,
        fastq1_input_files = TrimAdapters.trimmed_fastq1_files,
        fastq2_input_files = TrimAdapters.trimmed_fastq2_files,
        tar_star_reference = tar_star_reference
   }

    call Picard.RemoveDuplicatesFromBam as RemoveDuplicatesFromBam {
        input:
            input_ids = input_ids,
            aligned_bam_inputs = StarAlign.output_bam
    }

    call Picard.CollectMultipleMetricsMultiSample {
        input:
            aligned_bam_inputs = RemoveDuplicatesFromBam.output_bam,
            genome_ref_fasta = genome_ref_fasta,
            input_ids = input_ids
    }

    call CountAlignments.CountAlignments as CountAlignments {
        input:
            input_ids = input_ids,
            aligned_bam_inputs = RemoveDuplicatesFromBam.output_bam,
            annotation_gtf = annotations_gtf
    }

    call LoomUtils.SingleNucleusSmartSeq2LoomOutput as LoomOutput {
        input:
            input_ids = input_ids,
            input_names = input_names,
            pipeline_version = "SmartSeq2SingleNucleus_v~{pipeline_version}",
            input_id_metadata_field = input_id_metadata_field,
            input_name_metadata_field = input_name_metadata_field,
            alignment_summary_metrics = CollectMultipleMetricsMultiSample.alignment_summary_metrics,
            dedup_metrics = RemoveDuplicatesFromBam.dedup_metrics,
            gc_bias_summary_metrics = CollectMultipleMetricsMultiSample.gc_bias_summary_metrics,
            introns_counts = CountAlignments.intron_counts_out,
            exons_counts = CountAlignments.exon_counts_out,
            annotation_introns_added_gtf = annotations_gtf
    }

  ### Aggregate the Loom Files Directly ###
  call LoomUtils.AggregateSmartSeq2Loom as AggregateLoom {
    input:
      loom_input = LoomOutput.loom_output,
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
    Array[File] exon_intron_count_files = LoomOutput.exon_intron_counts
    Array[File] bam_files = RemoveDuplicatesFromBam.output_bam
    String pipeline_version_out = pipeline_version
  }
}
