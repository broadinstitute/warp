version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data processing according to the GATK Best Practices (June 2016)
## for human whole-genome and exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../../../../../tasks/broad/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../../../../tasks/broad/AggregatedBamQC.wdl" as AggregatedQC
import "../../../../../../tasks/broad/Utilities.wdl" as Utils

# WORKFLOW DEFINITION
workflow VUMCUnmappedBamToAlignedCram {
  input {
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    DragmapReference? dragmap_reference

    String cross_check_fingerprints_by = "READGROUP"
    Float lod_threshold = -20.0
    File haplotype_database_file
    Int preemptible_tries
    Int agg_preemptible_tries

    Float cutoff_for_large_rg_in_gb = 10.0
    Int reads_per_file = 48000000

    Boolean check_contaminant = true
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
    Boolean bin_base_qualities = true
    Boolean somatic = false
    Boolean use_bwa_mem = true
    Boolean perform_bqsr = true
    Boolean allow_empty_ref_alt = true

    Array[File] flowcell_unmapped_bams
    String sample_name

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    File calling_interval_list

    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File ref_alt
    File ref_sa
    File ref_amb
    File ref_bwt
    File ref_ann
    File ref_pac
    File? ref_str

    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices

    File dbsnp_vcf
    File dbsnp_vcf_index

    File evaluation_interval_list
  }

  SampleAndUnmappedBams sample_and_unmapped_bams = object {
    base_file_name: sample_name,
    final_gvcf_base_name: sample_name,
    flowcell_unmapped_bams: flowcell_unmapped_bams,
    sample_name: sample_name,
    unmapped_bam_suffix: ".bam"
  }

  ReferenceFasta reference_fasta = object {
    ref_dict: ref_dict,
    ref_fasta: ref_fasta,
    ref_fasta_index: ref_fasta_index,
    ref_alt: ref_alt,
    ref_sa: ref_sa,
    ref_amb: ref_amb,
    ref_bwt: ref_bwt,
    ref_ann: ref_ann,
    ref_pac: ref_pac,
    ref_str: ref_str
  }

  DNASeqSingleSampleReferences references = object {
    contamination_sites_ud: contamination_sites_ud,
    contamination_sites_bed: contamination_sites_bed,
    contamination_sites_mu: contamination_sites_mu,
    calling_interval_list: calling_interval_list,

    reference_fasta: reference_fasta,
    
    known_indels_sites_vcfs: known_indels_sites_vcfs,
    known_indels_sites_indices: known_indels_sites_indices,

    dbsnp_vcf: dbsnp_vcf,
    dbsnp_vcf_index: dbsnp_vcf_index,

    evaluation_interval_list: evaluation_interval_list,

    haplotype_database_file: haplotype_database_file
  }

  PapiSettings papi_settings = object {
    preemptible_tries: preemptible_tries,
    agg_preemptible_tries: agg_preemptible_tries
  }

  String recalibrated_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicates_marked.recalibrated"

  call ToBam.UnmappedBamToAlignedBam {
    input:
      sample_and_unmapped_bams    = sample_and_unmapped_bams,
      references                  = references,
      dragmap_reference           = dragmap_reference,
      papi_settings               = papi_settings,

      contamination_sites_ud = references.contamination_sites_ud,
      contamination_sites_bed = references.contamination_sites_bed,
      contamination_sites_mu = references.contamination_sites_mu,

      cross_check_fingerprints_by = cross_check_fingerprints_by,

      #we don't need to do CrossCheckFingerprints
      #haplotype_database_file     = references.haplotype_database_file,
      
      lod_threshold               = lod_threshold,
      recalibrated_bam_basename   = recalibrated_bam_basename,
      perform_bqsr                = perform_bqsr,
      use_bwa_mem                 = use_bwa_mem,
      unmap_contaminant_reads     = unmap_contaminant_reads,
      allow_empty_ref_alt         = allow_empty_ref_alt
  }

  call AggregatedQC.AggregatedBamQC {
    input:
      base_recalibrated_bam = UnmappedBamToAlignedBam.output_bam,
      base_recalibrated_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      base_name = sample_and_unmapped_bams.base_file_name,
      sample_name = sample_and_unmapped_bams.sample_name,
      recalibrated_bam_base_name = recalibrated_bam_basename,
      haplotype_database_file = references.haplotype_database_file,
      references = references,
      papi_settings = papi_settings
  }

  call Utils.ConvertToCram as ConvertToCram {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = sample_name,
      preemptible_tries = agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = UnmappedBamToAlignedBam.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = AggregatedBamQC.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = AggregatedBamQC.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = AggregatedBamQC.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = AggregatedBamQC.read_group_gc_bias_summary_metrics

    File selfSM = UnmappedBamToAlignedBam.selfSM
    Float contamination = UnmappedBamToAlignedBam.contamination

    File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = AggregatedBamQC.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = AggregatedBamQC.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = AggregatedBamQC.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = AggregatedBamQC.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = AggregatedBamQC.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
    File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    File? output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
    File output_cram_md5 = ConvertToCram.output_cram_md5
  }
  meta {
    allowNestedInputs: true
  }
}
