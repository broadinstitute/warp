version 1.0

import "../../dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl" as WholeGenomeGermlineSingleSample
import "../cram_to_unmapped_bams/CramToUnmappedBams.wdl" as ToUbams
import "../../../../structs/dna_seq/DNASeqStructs.wdl"

workflow WholeGenomeReprocessing {

  String pipeline_version = "3.3.5"

  input {
    File? input_cram
    File? input_bam
    File? output_map

    String sample_name
    String base_file_name
    String final_gvcf_base_name
    String unmapped_bam_suffix

    File? cram_ref_fasta
    File? cram_ref_fasta_index

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list

    String cloud_provider
  }

  call ToUbams.CramToUnmappedBams {
    input:
      input_cram = input_cram,
      input_bam = input_bam,
      ref_fasta = select_first([cram_ref_fasta, references.reference_fasta.ref_fasta]),
      ref_fasta_index = select_first([cram_ref_fasta_index, references.reference_fasta.ref_fasta_index]),
      base_file_name = base_file_name,
      output_map = output_map
  }

  SampleAndUnmappedBams sample_and_unmapped_bams = object {
     sample_name: sample_name,
     base_file_name: base_file_name,
     flowcell_unmapped_bams: CramToUnmappedBams.unmapped_bams,
     final_gvcf_base_name: final_gvcf_base_name,
     unmapped_bam_suffix: unmapped_bam_suffix
  }

  call WholeGenomeGermlineSingleSample.WholeGenomeGermlineSingleSample {
    input:
      sample_and_unmapped_bams = sample_and_unmapped_bams,
      references = references,
      scatter_settings = scatter_settings,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      cloud_provider = cloud_provider
  }

  output {
    Array[File] validation_report = CramToUnmappedBams.validation_report
    Array[File] unmapped_bams = CramToUnmappedBams.unmapped_bams

    Array[File] quality_yield_metrics = WholeGenomeGermlineSingleSample.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = WholeGenomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = WholeGenomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = WholeGenomeGermlineSingleSample.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = WholeGenomeGermlineSingleSample.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = WholeGenomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = WholeGenomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = WholeGenomeGermlineSingleSample.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = WholeGenomeGermlineSingleSample.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = WholeGenomeGermlineSingleSample.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = WholeGenomeGermlineSingleSample.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = WholeGenomeGermlineSingleSample.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = WholeGenomeGermlineSingleSample.read_group_gc_bias_summary_metrics

    File? cross_check_fingerprints_metrics = WholeGenomeGermlineSingleSample.cross_check_fingerprints_metrics

    File selfSM = WholeGenomeGermlineSingleSample.selfSM
    Float contamination = WholeGenomeGermlineSingleSample.contamination

    File calculate_read_group_checksum_md5 = WholeGenomeGermlineSingleSample.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = WholeGenomeGermlineSingleSample.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = WholeGenomeGermlineSingleSample.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = WholeGenomeGermlineSingleSample.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = WholeGenomeGermlineSingleSample.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = WholeGenomeGermlineSingleSample.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = WholeGenomeGermlineSingleSample.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = WholeGenomeGermlineSingleSample.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = WholeGenomeGermlineSingleSample.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = WholeGenomeGermlineSingleSample.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = WholeGenomeGermlineSingleSample.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = WholeGenomeGermlineSingleSample.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = WholeGenomeGermlineSingleSample.agg_quality_distribution_metrics

    File? fingerprint_summary_metrics = WholeGenomeGermlineSingleSample.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = WholeGenomeGermlineSingleSample.fingerprint_detail_metrics

    File wgs_metrics = WholeGenomeGermlineSingleSample.wgs_metrics
    File raw_wgs_metrics = WholeGenomeGermlineSingleSample.raw_wgs_metrics

    File duplicate_metrics = WholeGenomeGermlineSingleSample.duplicate_metrics
    File? output_bqsr_reports = WholeGenomeGermlineSingleSample.output_bqsr_reports

    File gvcf_summary_metrics = WholeGenomeGermlineSingleSample.gvcf_summary_metrics
    File gvcf_detail_metrics = WholeGenomeGermlineSingleSample.gvcf_detail_metrics

    File output_cram = WholeGenomeGermlineSingleSample.output_cram
    File output_cram_index = WholeGenomeGermlineSingleSample.output_cram_index
    File output_cram_md5 = WholeGenomeGermlineSingleSample.output_cram_md5

    File validate_cram_file_report = WholeGenomeGermlineSingleSample.validate_cram_file_report

    File output_vcf = WholeGenomeGermlineSingleSample.output_vcf
    File output_vcf_index = WholeGenomeGermlineSingleSample.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
