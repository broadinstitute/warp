version 1.0

import "JukeboxTasks.wdl" as Tasks
import "Qc.wdl" as QC

workflow QC {
  input {
    File agg_bam
    File agg_bam_index
    String base_file_name
    String base_file_name_sub

    Float agg_bam_size
    Float ref_size
    Int? additional_metrics_disk
    Float secure_disk_size_threshold

    References references
    ContaminationSites contamination_sites
    File wgs_coverage_interval_list

    File? picard_jar_override
    String gitc_path

    File monitoring_script
    Int VCF_disk_size
    Int additional_disk

    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
    String flow_order
  }
  Float dynamic_check_contamination_disk_size = agg_bam_size + ref_size + additional_metrics_disk
  Float check_contamination_disk_size = if dynamic_check_contamination_disk_size > secure_disk_size_threshold then dynamic_check_contamination_disk_size else secure_disk_size_threshold

  # Preprocessing step for VerifyBamID, where we realign the reads to the ref and alt haplotypes and extract a bamout.
  # --alleles ~{contamination_sites_vcf} activates the "genotype-given-allele" mode. This way, even in the case where the alt alleles have been pushed out
  # and the pileup is unable to see them, the alt haplotype is incorporated into the de bruijn graph and so the read likelihood engine should be able to figure out
  # that the read aligns just as well to the alt haplotype.
  #
  # We switched from the pairHMM read likelihood to engine to the flow-based read likelihood engine
  # in November 2021.

  # This vcf contains the SNPs sites with high allele frequency---the same sites as contamination_sites_bed.
  # We need to convert the bed to a vcf because we need alleles, not just the contig and positions, for the genotype-given-allele mode (i.e. --allele argument).

  String hc_contamination_extra_args = "--bam-writer-type NO_HAPLOTYPES --alleles " + contamination_sites.contamination_sites_vcf
  call Tasks.HaplotypeCaller as HaplotypeCallerForContamination {
    input:
      input_bam_list = [agg_bam],
      interval_list = contamination_sites.contamination_sites_vcf,
      vcf_basename = base_file_name,
      references = references,
      disk_size = ceil((agg_bam_size + VCF_disk_size) + ref_size + additional_disk),
      gitc_path = gitc_path,
      extra_args = hc_contamination_extra_args,
      make_gvcf = false,
      memory_gb = 12,
      make_bamout = true
  }

  File contamination_sites_ud = contamination_sites.contamination_sites_path + flow_order + ".UD"
  File contamination_sites_bed = contamination_sites.contamination_sites_path + flow_order + ".bed"
  File contamination_sites_mu = contamination_sites.contamination_sites_path + flow_order + ".mu"

  # Estimate level of cross-sample contamination
  call Tasks.CheckContamination {
    input:
      input_bam = select_first([HaplotypeCallerForContamination.bamout]),
      input_bam_index = select_first([HaplotypeCallerForContamination.bamout_index]),
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      references = references,
      output_prefix = base_file_name,
      disk_size = ceil(check_contamination_disk_size)
  }

  Float dynamic_statistics_disk_size = agg_bam_size + ref_size + ( additional_metrics_disk * 2 )
  Int statistics_disk_size = if dynamic_statistics_disk_size > secure_disk_size_threshold then ceil(dynamic_statistics_disk_size) else ceil(secure_disk_size_threshold)

  call Tasks.CollectDuplicateMetrics {
    input:
      input_bam = agg_bam,
      metrics_filename = base_file_name_sub + ".duplicate_metrics",
      disk_size_gb = statistics_disk_size,
      gitc_path = gitc_path,
      jar_override = picard_jar_override,
  }

  call QC.CollectQualityYieldMetrics {
    input:
      input_bam = agg_bam,
      metrics_filename = base_file_name_sub + ".unmapped.quality_yield_metrics",
  }

  # QC the sample WGS metrics (stringent thresholds)
  call Tasks.CollectWgsMetrics {
    input:
      input_bam = agg_bam,
      input_bam_index = agg_bam_index,
      metrics_filename = base_file_name_sub + ".wgs_metrics",
      references = references,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      disk_size = statistics_disk_size,
      gitc_path = gitc_path,
      jar_override = picard_jar_override
  }

  Int default_raw_wgs_memory_size = 12
  Int increased_raw_wgs_memory_size = 30
  Int input_bam_size_threshold = 200
  Int raw_wgs_memory_size = if ceil(agg_bam_size) > input_bam_size_threshold then increased_raw_wgs_memory_size else default_raw_wgs_memory_size

  # QC the sample raw WGS metrics (common thresholds)
  call Tasks.CollectRawWgsMetrics {
    input:
      input_bam = agg_bam,
      input_bam_index = agg_bam_index,
      metrics_filename = base_file_name_sub + ".raw_wgs_metrics",
      references = references,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      disk_size = statistics_disk_size,
      memory_size = raw_wgs_memory_size,
      gitc_path = gitc_path,
      jar_override = picard_jar_override
  }

  # QC the final BAM some more (no such thing as too much QC)
  call Tasks.CollectAggregationMetrics {
    input:
      input_bam = agg_bam,
      input_bam_index = agg_bam_index,
      output_bam_prefix = base_file_name_sub,
      references = references,
      disk_size = statistics_disk_size,
      gitc_path = gitc_path,
      jar_override = picard_jar_override
  }

  # Check whether the data has massively high duplication or chimerism rates

  File dup_metrics = CollectDuplicateMetrics.duplicate_metrics
  call QC.CheckPreValidation {
    input:
      duplication_metrics = dup_metrics,
      chimerism_metrics = CollectAggregationMetrics.alignment_summary_metrics,
      max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample,
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