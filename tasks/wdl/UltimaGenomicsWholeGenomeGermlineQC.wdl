version 1.0

import "../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineTasks.wdl" as Tasks
import "../../tasks/wdl/Qc.wdl" as QC

workflow UltimaGenomicsWholeGenomeGermlineQC {
  input {
    File agg_bam
    File agg_bam_index
    String base_file_name
    String base_file_name_sub

    References references
    ContaminationSites contamination_sites
    File wgs_coverage_interval_list

    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
    String flow_order
  }

  # Preprocessing step for VerifyBamID, where we realign the reads to the ref and alt haplotypes and extract a bamout.
  # --alleles ~{contamination_sites_vcf} activates the "genotype-given-allele" mode. This way, even in the case where the alt alleles have been pushed out
  # and the pileup is unable to see them, the alt haplotype is incorporated into the de bruijn graph and so the read likelihood engine should be able to figure out
  # that the read aligns just as well to the alt haplotype.
  #
  # We switched from the pairHMM read likelihood to engine to the flow-based read likelihood engine
  # in November 2021.

  # This vcf contains the SNPs sites with high allele frequency---the same sites as contamination_sites_bed.
  # We need to convert the bed to a vcf because we need alleles, not just the contig and positions, for the genotype-given-allele mode (i.e. --allele argument).

  String hc_contamination_extra_args = "--alleles " + contamination_sites.contamination_sites_vcf
  call Tasks.HaplotypeCaller as HaplotypeCallerForContamination {
    input:
      input_bam_list                = [agg_bam],
      input_bam_index_list          = [agg_bam_index],
      interval_list                 = contamination_sites.contamination_sites_vcf,
      vcf_basename                  = base_file_name,
      references                    = references,
      contamination_extra_args      = hc_contamination_extra_args,
      make_bamout                   = true
  }

  String contamination_sites_ud   = contamination_sites.contamination_sites_path + flow_order + ".UD"
  String contamination_sites_bed  = contamination_sites.contamination_sites_path + flow_order + ".bed"
  String contamination_sites_mu   = contamination_sites.contamination_sites_path + flow_order + ".mu"

  # Estimate level of cross-sample contamination
  call Tasks.CheckContamination {
    input:
      input_bam               = HaplotypeCallerForContamination.bamout,
      input_bam_index         = HaplotypeCallerForContamination.bamout_index,
      contamination_sites_ud  = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu  = contamination_sites_mu,
      references              = references,
      output_prefix           = base_file_name
  }

  call Tasks.CollectDuplicateMetrics {
    input:
      input_bam        = agg_bam,
      metrics_filename = base_file_name_sub + ".duplicate_metrics",
      references       = references
  }

  call QC.CollectQualityYieldMetrics {
    input:
      input_bam        = agg_bam,
      metrics_filename = base_file_name_sub + ".unmapped.quality_yield_metrics",
  }

  # QC the sample WGS metrics (stringent thresholds)
  call Tasks.CollectWgsMetrics {
    input:
      input_bam                   = agg_bam,
      input_bam_index             = agg_bam_index,
      metrics_filename            = base_file_name_sub + ".wgs_metrics",
      references                  = references,
      wgs_coverage_interval_list  = wgs_coverage_interval_list
  }

  # QC the sample raw WGS metrics (common thresholds)
  call Tasks.CollectRawWgsMetrics {
    input:
      input_bam                   = agg_bam,
      input_bam_index             = agg_bam_index,
      metrics_filename            = base_file_name_sub + ".raw_wgs_metrics",
      references                  = references,
      wgs_coverage_interval_list  = wgs_coverage_interval_list
  }

  # QC the final BAM some more (no such thing as too much QC)
  call Tasks.CollectAggregationMetrics {
    input:
      input_bam           = agg_bam,
      input_bam_index     = agg_bam_index,
      output_bam_prefix   = base_file_name_sub,
      references          = references
  }

  # Check whether the data has massively high duplication or chimerism rates

  File dup_metrics = CollectDuplicateMetrics.duplicate_metrics
  call QC.CheckPreValidation {
    input:
      duplication_metrics                  = dup_metrics,
      chimerism_metrics                    = CollectAggregationMetrics.alignment_summary_metrics,
      max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample   = max_chimerism_in_reasonable_sample,
  }

  output {
    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination
    File quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics
    File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics
    File duplicate_metrics = CollectDuplicateMetrics.duplicate_metrics
    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
    File? agg_alignment_summary_pdf = CollectAggregationMetrics.alignment_summary_pdf
    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics
    Float duplication_rate = CheckPreValidation.duplication_rate
    Float chimerism_rate = CheckPreValidation.chimerism_rate
    Boolean is_outlier_data = CheckPreValidation.is_outlier_data
  }
}