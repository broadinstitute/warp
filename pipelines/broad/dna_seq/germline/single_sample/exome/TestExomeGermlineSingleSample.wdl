version 1.0

import "../../../../../..pipelines/broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as ExomeGermlineSingleSample
import "../../../../../../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample

workflow TestExomeGermlineSingleSample {

  input {
    PapiSettings papi_settings
    SampleAndUnmappedBams sample_and_unmapped_bams
    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index
    File target_interval_list
    File bait_interval_list
    String bait_set_name
    Boolean provide_bam_output = false

    String truth_cloud_path
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
  }

  # Get a list of metrics files to test from the truth bucket


  Array[String] truth_metrics =  prefix(truth_cloud_path, metrics_files_to_test)
  File truth_cram = truth_cloud_path + basename(ExomeGermlineSingleSample.output_cram)
  File truth_crai = truth_cloud_path + basename(ExomeGermlineSingleSample.output_cram_index)
  File truth_gvcf = truth_cloud_path + basename(ExomeGermlineSingleSample.output_vcf)

  call VerifyGermlineSingleSample.VerifyGermlineSingleSample {
    input:
    #truth_metrics = # Array[File] truth_metrics
    truth_cram = truth_cram,                                       # File truth_cram
    truth_crai = truth_crai,                                       # File truth_crai
    truth_gvcf = truth_gvcf,                                       # File truth_gvcf
    #test_metrics = # Array[File] test_metric
    test_cram = ExomeGermlineSingleSample.output_cram,             # File test_cram
    test_crai = ExomeGermlineSingleSample.output_cram_index,       # File test_crai
    test_gvcf = ExomeGermlineSingleSample.output_vcf              # File test_gvcf
  }

  # Outputs are intentionally  left out to automatically propagate all the outputs all tasks (wihout needing to inddividually define the  \outputs)


 #output {
 #  Array[File] quality_yield_metrics = UnmappedBamToAlignedBam.quality_yield_metrics

 #  Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_pdf
 #  Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_metrics
 #  Array[File] unsorted_read_group_insert_size_histogram_pdf = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_histogram_pdf
 #  Array[File] unsorted_read_group_insert_size_metrics = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_metrics
 #  Array[File] unsorted_read_group_quality_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_pdf
 #  Array[File] unsorted_read_group_quality_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_metrics
 #  Array[File] unsorted_read_group_quality_distribution_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_pdf
 #  Array[File] unsorted_read_group_quality_distribution_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_metrics

 #  File read_group_alignment_summary_metrics = AggregatedBamQC.read_group_alignment_summary_metrics

 #  File? cross_check_fingerprints_metrics = UnmappedBamToAlignedBam.cross_check_fingerprints_metrics

 #  File selfSM = UnmappedBamToAlignedBam.selfSM
 #  Float contamination = UnmappedBamToAlignedBam.contamination

 #  File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5

 #  File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
 #  File agg_bait_bias_detail_metrics = AggregatedBamQC.agg_bait_bias_detail_metrics
 #  File agg_bait_bias_summary_metrics = AggregatedBamQC.agg_bait_bias_summary_metrics
 #  File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
 #  File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
 #  File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
 #  File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
 #  File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
 #  File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
 #  File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

 #  File? fingerprint_summary_metrics = AggregatedBamQC.fingerprint_summary_metrics
 #  File? fingerprint_detail_metrics = AggregatedBamQC.fingerprint_detail_metrics

 #  File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
 #  File output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

 #  File gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
 #  File gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics

 #  File hybrid_selection_metrics = CollectHsMetrics.metrics

 #  File? output_bam = provided_output_bam
 #  File? output_bam_index = provided_output_bam_index

 #  File output_cram = BamToCram.output_cram
 #  File output_cram_index = BamToCram.output_cram_index
 #  File output_cram_md5 = BamToCram.output_cram_md5

 #  File validate_cram_file_report = BamToCram.validate_cram_file_report

 #  File output_vcf = BamToGvcf.output_vcf
 #  File output_vcf_index = BamToGvcf.output_vcf_index
 #}
}

task get_metrics_files {
  input {
    String truth_cloud_path
  }

  command {
    for i in $(gsutil ls $truth_cloud_path); do
  }

  runtime {
  }
  output {
    Array[String] metrics_files_to_test =
  }
}