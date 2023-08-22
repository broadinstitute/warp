version 1.0

## Copyright Broad Institute/VUMC, 2018/2022
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in FASTQ format
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../../../../../tasks/vumc_biostatistics/PairedFastQsToUnmappedBAM.wdl" as ToUnmappedBam
import "../../../../../../pipelines/broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as BroadPipeline
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow VUMCExomeGermlineSingleSampleFromFastq {

  String pipeline_version = "3.1.10"

  input {
    # Optional for VUMC pipeline
    String sample_name 
    String fastq_1 
    String fastq_2 
    String readgroup_name 
    String library_name 
    String platform_unit 
    String platform_name 
    String? run_date 
    String? sequencing_center 

    # Optional for BROAD pipeline
    PapiSettings papi_settings
    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File target_interval_list
    File bait_interval_list
    String bait_set_name

    Boolean provide_bam_output = false
  }

  # Convert pair of FASTQs to uBAM
  call ToUnmappedBam.PairedFastQsToUnmappedBAM {
    input:
      sample_name = sample_name,
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,
      readgroup_name = readgroup_name,
      library_name = library_name,
      platform_unit = platform_unit,
      run_date = run_date,
      platform_name = platform_name,
      sequencing_center = sequencing_center,
  }

  SampleAndUnmappedBams sample_and_unmapped_bams = object {
    base_file_name: sample_name,
    final_gvcf_base_name: sample_name,
    flowcell_unmapped_bams: [ PairedFastQsToUnmappedBAM.output_unmapped_bam ],
    sample_name: sample_name,
    unmapped_bam_suffix: ".bam"
  }

  call BroadPipeline.ExomeGermlineSingleSample as broad {
    input:
      papi_settings = papi_settings,
      sample_and_unmapped_bams = sample_and_unmapped_bams,
      references = references,
      scatter_settings = scatter_settings,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      target_interval_list = target_interval_list,
      bait_interval_list = bait_interval_list,
      bait_set_name = bait_set_name,
      provide_bam_output = provide_bam_output
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = broad.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = broad.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = broad.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = broad.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = broad.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = broad.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = broad.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = broad.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = broad.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = broad.read_group_alignment_summary_metrics

    File? cross_check_fingerprints_metrics = broad.cross_check_fingerprints_metrics

    File selfSM = broad.selfSM
    Float contamination = broad.contamination

    File calculate_read_group_checksum_md5 = broad.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = broad.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = broad.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = broad.agg_bait_bias_summary_metrics
    File agg_insert_size_histogram_pdf = broad.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = broad.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = broad.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = broad.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = broad.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = broad.agg_quality_distribution_metrics
    File agg_error_summary_metrics = broad.agg_error_summary_metrics

    File? fingerprint_summary_metrics = broad.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = broad.fingerprint_detail_metrics

    File duplicate_metrics = broad.duplicate_metrics
    File? output_bqsr_reports = broad.output_bqsr_reports

    File gvcf_summary_metrics = broad.gvcf_summary_metrics
    File gvcf_detail_metrics = broad.gvcf_detail_metrics

    File hybrid_selection_metrics = broad.hybrid_selection_metrics

    File? output_bam = broad.output_bam
    File? output_bam_index = broad.output_bam_index

    File output_cram = broad.output_cram
    File output_cram_index = broad.output_cram_index
    File output_cram_md5 = broad.output_cram_md5

    File validate_cram_file_report = broad.validate_cram_file_report

    File output_vcf = broad.output_vcf
    File output_vcf_index = broad.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}

