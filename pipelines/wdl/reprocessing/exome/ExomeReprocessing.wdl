version 1.0

import "../../dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as ExomeGermlineSingleSample
import "../cram_to_unmapped_bams/CramToUnmappedBams.wdl" as ToUbams
import "../../../../structs/dna_seq/DNASeqStructs.wdl"

workflow ExomeReprocessing {


  String pipeline_version = "3.3.7"

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

    File target_interval_list
    File bait_interval_list
    String bait_set_name

    String cloud_provider
  }

  call ToUbams.CramToUnmappedBams {
    input:
      input_cram = input_cram,
      input_bam = input_bam,
      ref_fasta = select_first([cram_ref_fasta, references.reference_fasta.ref_fasta]),
      ref_fasta_index = select_first([cram_ref_fasta_index, references.reference_fasta.ref_fasta_index]),
      output_map = output_map,
      base_file_name = base_file_name,
      unmapped_bam_suffix = unmapped_bam_suffix
  }

  SampleAndUnmappedBams sample_and_unmapped_bams = object {
     sample_name: sample_name,
     base_file_name: base_file_name,
     flowcell_unmapped_bams: CramToUnmappedBams.unmapped_bams,
     final_gvcf_base_name: final_gvcf_base_name,
     unmapped_bam_suffix: unmapped_bam_suffix
  }

  call ExomeGermlineSingleSample.ExomeGermlineSingleSample {
    input:
      sample_and_unmapped_bams = sample_and_unmapped_bams,
      references = references,
      scatter_settings = scatter_settings,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings,
      target_interval_list = target_interval_list,
      bait_interval_list = bait_interval_list,
      bait_set_name = bait_set_name,
      cloud_provider = cloud_provider
  }

  output {
    Array[File] validation_report = CramToUnmappedBams.validation_report
    Array[File] unmapped_bams = CramToUnmappedBams.unmapped_bams

    Array[File] quality_yield_metrics = ExomeGermlineSingleSample.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = ExomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = ExomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = ExomeGermlineSingleSample.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = ExomeGermlineSingleSample.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = ExomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = ExomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = ExomeGermlineSingleSample.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = ExomeGermlineSingleSample.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = ExomeGermlineSingleSample.read_group_alignment_summary_metrics

    File? cross_check_fingerprints_metrics = ExomeGermlineSingleSample.cross_check_fingerprints_metrics

    File selfSM = ExomeGermlineSingleSample.selfSM
    Float contamination = ExomeGermlineSingleSample.contamination

    File calculate_read_group_checksum_md5 = ExomeGermlineSingleSample.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = ExomeGermlineSingleSample.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = ExomeGermlineSingleSample.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = ExomeGermlineSingleSample.agg_bait_bias_summary_metrics
    File agg_insert_size_histogram_pdf = ExomeGermlineSingleSample.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = ExomeGermlineSingleSample.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = ExomeGermlineSingleSample.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = ExomeGermlineSingleSample.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = ExomeGermlineSingleSample.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = ExomeGermlineSingleSample.agg_quality_distribution_metrics
    File agg_error_summary_metrics = ExomeGermlineSingleSample.agg_error_summary_metrics

    File? fingerprint_summary_metrics = ExomeGermlineSingleSample.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = ExomeGermlineSingleSample.fingerprint_detail_metrics

    File duplicate_metrics = ExomeGermlineSingleSample.duplicate_metrics
    File? output_bqsr_reports = ExomeGermlineSingleSample.output_bqsr_reports

    File gvcf_summary_metrics = ExomeGermlineSingleSample.gvcf_summary_metrics
    File gvcf_detail_metrics = ExomeGermlineSingleSample.gvcf_detail_metrics

    File hybrid_selection_metrics = ExomeGermlineSingleSample.hybrid_selection_metrics

    File output_cram = ExomeGermlineSingleSample.output_cram
    File output_cram_index = ExomeGermlineSingleSample.output_cram_index
    File output_cram_md5 = ExomeGermlineSingleSample.output_cram_md5

    File validate_cram_file_report = ExomeGermlineSingleSample.validate_cram_file_report

    File output_vcf = ExomeGermlineSingleSample.output_vcf
    File output_vcf_index = ExomeGermlineSingleSample.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
