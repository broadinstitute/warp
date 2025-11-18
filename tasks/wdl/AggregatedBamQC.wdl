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

import "../../tasks/wdl/Qc.wdl" as QC
import "../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow AggregatedBamQC {
input {
    File base_recalibrated_bam
    File base_recalibrated_bam_index
    String base_name
    String sample_name
    String recalibrated_bam_base_name
    File haplotype_database_file
    DNASeqSingleSampleReferences references
    PapiSettings papi_settings
    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index
  }

  # QC the final BAM (consolidated after scattered BQSR)
  call QC.CollectReadgroupBamQualityMetrics as CollectReadgroupBamQualityMetrics {
    input:
      input_bam = base_recalibrated_bam,
      input_bam_index = base_recalibrated_bam_index,
      output_bam_prefix = base_name + ".readgroup",
      ref_dict = references.reference_fasta.ref_dict,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # QC the final BAM some more (no such thing as too much QC)
  call QC.CollectAggregationMetrics as CollectAggregationMetrics {
    input:
      input_bam = base_recalibrated_bam,
      input_bam_index = base_recalibrated_bam_index,
      output_bam_prefix = base_name,
      ref_dict = references.reference_fasta.ref_dict,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  if (defined(haplotype_database_file) && defined(fingerprint_genotypes_file)) {
    # Check the sample BAM fingerprint against the sample array
    call QC.CheckFingerprintTask as CheckFingerprintTask {
      input:
        input_bam = base_recalibrated_bam,
        input_bam_index = base_recalibrated_bam_index,
        genotypes = select_first([fingerprint_genotypes_file]),
        genotypes_index = fingerprint_genotypes_index,
        expected_sample_alias = sample_name,
        output_basename = base_name,
        haplotype_database_file = haplotype_database_file,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Generate a checksum per readgroup in the final BAM
  call QC.CalculateReadGroupChecksum as CalculateReadGroupChecksum {
    input:
      input_bam = base_recalibrated_bam,
      input_bam_index = base_recalibrated_bam_index,
      read_group_md5_filename = recalibrated_bam_base_name + ".bam.read_group_md5",
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  output {
    File read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics
    File read_group_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf
    File read_group_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics

    File calculate_read_group_checksum_md5 = CalculateReadGroupChecksum.md5_file

    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
    File agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
    File agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
    File agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics
    File agg_error_summary_metrics = CollectAggregationMetrics.error_summary_metrics

    File? fingerprint_summary_metrics = CheckFingerprintTask.summary_metrics
    File? fingerprint_detail_metrics = CheckFingerprintTask.detail_metrics
  }
  meta {
    allowNestedInputs: true
  }
}
