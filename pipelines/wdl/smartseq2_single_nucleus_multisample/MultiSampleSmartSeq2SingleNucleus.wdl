version 1.0

import "../../../tasks/wdl/CheckInputs.wdl" as CheckInputs
import "../../../tasks/wdl/TrimAdapters.wdl" as TrimAdapters
import "../../../tasks/wdl/StarAlign.wdl" as StarAlign
import "../../../tasks/wdl/Picard.wdl" as Picard
import "../../../tasks/wdl/FeatureCounts.wdl" as CountAlignments
import "../../../tasks/wdl/H5adUtils.wdl" as H5adUtils
import "../../../tasks/wdl/Utilities.wdl" as utils

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

      String cloud_provider
  }

  String ubuntu_docker = "ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
  String gcp_ubuntu_docker_prefix = "gcr.io/gcp-runtimes/"
  String acr_ubuntu_docker_prefix = "dsppipelinedev.azurecr.io/"
  String ubuntu_docker_prefix = if cloud_provider == "gcp" then gcp_ubuntu_docker_prefix else acr_ubuntu_docker_prefix

  # make sure either gcp or azr is supplied as cloud_provider input
  if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
      call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
          input:
              message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
      }
  }

  # Version of this pipeline
  String pipeline_version = "2.2.3"

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
  
  call StarAlign.STARGenomeRefVersion as ReferenceCheck {
    input:
      tar_star_reference = tar_star_reference,
      ubuntu_docker_path = ubuntu_docker_prefix + ubuntu_docker
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

    call H5adUtils.SingleNucleusSmartSeq2H5adOutput as H5adOutput {
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

  ### Aggregate the H5ad Files Directly ###
  call H5adUtils.AggregateSmartSeq2H5ad as AggregateH5ad {
    input:
      h5ad_input = H5adOutput.h5ad_output,
      pipeline_version = pipeline_version,
      batch_id = batch_id
  }



  ### Pipeline output ###
  output {
    # h5ad output, exon/intron count tsv files and the aligned bam files
    File h5ad_output = AggregateH5ad.h5ad_output_file
    File genomic_reference_version = ReferenceCheck.genomic_ref_version
    Array[File] exon_intron_count_files = H5adOutput.exon_intron_counts
    Array[File] bam_files = RemoveDuplicatesFromBam.output_bam
    String pipeline_version_out = pipeline_version
  }
}
