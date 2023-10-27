version 1.0

import "../../../tasks/skylab/CheckInputs.wdl" as CheckInputs
import "../../../tasks/skylab/TrimAdapters.wdl" as TrimAdapters
import "../../../tasks/skylab/StarAlign.wdl" as StarAlign
import "../../../tasks/skylab/Picard.wdl" as Picard
import "../../../tasks/skylab/FeatureCounts.wdl" as CountAlignments
import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils
import "../../../tasks/broad/Utilities.wdl" as Utils


workflow MultiSampleSmartSeq2SingleNucleus {
  meta {
    description: "The MultiSampleSmartSeq2SingleNucleus pipeline runs multiple snSS2 samples in a single pipeline invocation"
    allowNestedInputs: true
  }

  input {
      # Cloud provider
      String cloud_provider = "gcp"

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
  String pipeline_version = "1.2.25"

  # Docker images
  String picard_cloud_docker = "picard-cloud:2.26.10"
  String alpine_docker = "alpine-bash:latest"
  String ubuntu_docker = "ubuntu_16_0_4:latest"
  String ea_utils_docker = "ea-utils:1.0.0-1.04.807-1659990665"
  String star_docker = "star:1.0.0-2.7.9a-1658781884"
  String subread_docker = "subread:1.0.0-2.0.1-1689097353"
  String pytools_docker = "pytools:1.0.0-1661263730"

  #TODO how do we handle these?
  String gcp_alpine_docker_prefix = "bashell/"
  String acr_alpine_docker_prefix = "dsppipelinedev.azurecr.io/"
  String alpine_docker_prefix = if cloud_provider == "gcp" then gcp_alpine_docker_prefix else acr_alpine_docker_prefix

  String ubuntu_docker = "ubuntu_16_0_4:latest"
  String gcp_ubuntu_docker_prefix = "gcr.io/gcp-runtimes/"
  String acr_ubuntu_docker_prefix = "dsppipelinedev.azurecr.io/"
  String ubuntu_docker_prefix = if cloud_provider == "gcp" then gcp_ubuntu_docker_prefix else acr_ubuntu_docker_prefix


  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String acr_docker_prefix = "dsppipelinedev.azurecr.io/"

  # choose docker prefix based on cloud provider
  String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

  # make sure either gcp or azr is supplied as cloud_provider input
  if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
      call Utils.ErrorWithMessage as ErrorMessageIncorrectInput {
          input:
              message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
      }
  }

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
         paired_end = true,
         alpine_docker_path = alpine_docker_prefix + alpine_docker

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
         adapter_list = adapter_list,
         ea_utils_docker_path = docker_prefix + ea_utils_docker
   }

   call StarAlign.StarAlignFastqMultisample as StarAlign {
      input:
        input_ids = input_ids,
        fastq1_input_files = TrimAdapters.trimmed_fastq1_files,
        fastq2_input_files = TrimAdapters.trimmed_fastq2_files,
        tar_star_reference = tar_star_reference,
        star_docker_path = docker_prefix + star_docker

   }

    call Picard.RemoveDuplicatesFromBam as RemoveDuplicatesFromBam {
        input:
            input_ids = input_ids,
            aligned_bam_inputs = StarAlign.output_bam,
            picard_docker_path = docker_prefix + picard_cloud_docker
    }

    call Picard.CollectMultipleMetricsMultiSample {
        input:
            aligned_bam_inputs = RemoveDuplicatesFromBam.output_bam,
            genome_ref_fasta = genome_ref_fasta,
            input_ids = input_ids,
            picard_docker_path = docker_prefix + picard_cloud_docker
    }

    call CountAlignments.CountAlignments as CountAlignments {
        input:
            input_ids = input_ids,
            aligned_bam_inputs = RemoveDuplicatesFromBam.output_bam,
            annotation_gtf = annotations_gtf,
            subread_docker_path = docker_prefix + subread_docker
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
            annotation_introns_added_gtf = annotations_gtf,
            pytools_docker_path = docker_prefix + pytools_docker
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
      pipeline_version = "MultiSampleSmartSeq2SingleNucleus_v~{pipeline_version}",
      pytools_docker_path = docker_prefix + pytools_docker
  }



  ### Pipeline output ###
  output {
    # loom output, exon/intron count tsv files and the aligned bam files
    File loom_output = AggregateLoom.loom_output_file
    File genomic_reference_version = ReferenceCheck.genomic_ref_version
    Array[File] exon_intron_count_files = LoomOutput.exon_intron_counts
    Array[File] bam_files = RemoveDuplicatesFromBam.output_bam
    String pipeline_version_out = pipeline_version
  }
}
