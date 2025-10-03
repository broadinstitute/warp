version 1.0

import "../../../../../../tasks/wdl/JointGenotypingTasks.wdl" as Tasks
import "https://raw.githubusercontent.com/broadinstitute/gatk/4.5.0.0/scripts/vcf_site_level_filtering_wdl/JointVcfFiltering.wdl" as Filtering
import "../../../../../../tasks/wdl/UltimaGenomicsGermlineFilteringThreshold.wdl" as FilteringThreshold


# Joint Genotyping for hg38 Whole Genomes (has not been tested on hg19) sequenced with Ultima sequencer.
# Joint genotyping performed with GenomicsDB, Filtering performed by GATK adaptable site level filtering pipeline.
# Filtering thershold chosen with Ultima's variant calling package by maximizing F1 for a known truth sample.
# For choosing a filtering threshold (where on the ROC curve to filter) a sample with truth data is required.
workflow UltimaGenomicsJointGenotyping {

  String pipeline_version = "1.2.3"

  input {
    File unpadded_intervals_file

    String callset_name
    File sample_name_map

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    Int small_disk = 100
    Int medium_disk = 200
    Int large_disk = 1000
    Int huge_disk = 2000

    File haplotype_database

    File eval_interval_list

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69

    #inputs for threshold filtering
    String truth_sample_name
    File truth_vcf
    File truth_vcf_index
    File truth_highconf_intervals
    String call_sample_name
    File ref_fasta_sdf
    File runs_file
    Array[File] annotation_intervals
    String flow_order

    #inputs for training and applying filter model
    Array[String] snp_annotations
    Array[String] indel_annotations
    String model_backend
    String snp_resource_args = "--resource:hapmap,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz --resource:omni,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz --resource:1000G,training=true,calibration=false gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    String indel_resource_args = "--resource:mills,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

    Int? top_level_scatter_count
    Boolean? gather_vcfs
    Float unbounded_scatter_count_scale_factor = 0.15
    Boolean cross_check_fingerprints = true
    Boolean scatter_cross_check_fingerprints = false
  }

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000

  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  # We allow overriding this default behavior for testing / special requests.
  Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])

  Int unbounded_scatter_count = select_first([top_level_scatter_count, round(unbounded_scatter_count_scale_factor * num_gvcfs)])
  Int scatter_count = if unbounded_scatter_count > 2 then unbounded_scatter_count else 2 #I think weird things happen if scatterCount is 1 -- IntervalListTools is noop?

  call Tasks.CheckSamplesUnique {
    input:
      sample_name_map = sample_name_map
  }

  call Tasks.SplitIntervalList {
    input:
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size_gb = small_disk,
      sample_names_unique_done = CheckSamplesUnique.samples_unique
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call Tasks.ImportGVCFs {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        workspace_dir_name = "genomicsdb",
        disk_size_gb = medium_disk,
        batch_size = 50
    }

    call Tasks.GenotypeGVCFs {
        input:
          workspace_tar = ImportGVCFs.output_genomicsdb,
          interval = unpadded_intervals[idx],
          output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf,
          disk_size_gb = medium_disk,
          keep_combined_raw_annotations = true,
          additional_annotation = "RawGtCount"
    }

    #TODO: if this task becomes expensive or slow we should combine the functionality into GenotypeGVCFs in GATK
    call Tasks.CalculateAverageAnnotations {
      input:
        vcf = GenotypeGVCFs.output_vcf,
        vcf_index = GenotypeGVCFs.output_vcf_index
    }

    call Tasks.HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = CalculateAverageAnnotations.output_vcf,
        vcf_index = CalculateAverageAnnotations.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        disk_size_gb = medium_disk
    }
  }

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size_gb = medium_disk
  }

  call Filtering.JointVcfFiltering as TrainAndApplyFilteringModelSNPs {
    input:
      input_vcfs = CalculateAverageAnnotations.output_vcf,
      input_vcf_idxs = CalculateAverageAnnotations.output_vcf_index,
      sites_only_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_vcf_idx = SitesOnlyGatherVcf.output_vcf_index,
      annotations = snp_annotations,
      resource_args = snp_resource_args,
      model_backend = model_backend,
      output_prefix = callset_name,
      extract_extra_args = "--mode SNP",
      train_extra_args = "--mode SNP",
      score_extra_args = "--mode SNP",
      gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }

  call Filtering.JointVcfFiltering as TrainAndApplyFilteringModelINDELs {
    input:
      input_vcfs = TrainAndApplyFilteringModelSNPs.scored_vcfs,
      input_vcf_idxs = TrainAndApplyFilteringModelSNPs.scored_vcf_idxs,
      sites_only_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_vcf_idx = SitesOnlyGatherVcf.output_vcf_index,
      annotations = indel_annotations,
      resource_args = indel_resource_args,
      model_backend = model_backend,
      output_prefix = callset_name,
      extract_extra_args = "--mode INDEL",
      train_extra_args = "--mode INDEL",
      score_extra_args = "--mode INDEL",
      gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }

  call FilteringThreshold.ExtractOptimizeSingleSample as FindFilteringThresholdAndFilter {
    input:
      input_vcf = TrainAndApplyFilteringModelINDELs.scored_vcfs,
      input_vcf_index = TrainAndApplyFilteringModelINDELs.scored_vcf_idxs,
      base_file_name = callset_name,
      call_sample_name = call_sample_name,
      truth_vcf = truth_vcf,
      truth_vcf_index = truth_vcf_index,
      truth_highconf_intervals = truth_highconf_intervals,
      truth_sample_name = truth_sample_name,
      flow_order = flow_order,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_fasta_sdf = ref_fasta_sdf,
      runs_file = runs_file,
      annotation_intervals = annotation_intervals,
      medium_disk = medium_disk
  }

  scatter (idx in range(length(TrainAndApplyFilteringModelINDELs.scored_vcfs))) {
    # For large callsets we need to collect metrics from the shards and gather them later.
    if (!is_small_callset) {
      call Tasks.CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = FindFilteringThresholdAndFilter.output_vcf[idx],
          input_vcf_index = FindFilteringThresholdAndFilter.output_vcf_index[idx],
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size_gb = medium_disk
      }
    }
  }

  # For small callsets we can gather the VCF shards and then collect metrics on it.
  if (is_small_callset) {
    call Tasks.GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs = FindFilteringThresholdAndFilter.output_vcf,
        output_vcf_name = callset_name + ".vcf.gz",
        disk_size_gb = huge_disk
    }

    call Tasks.CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size_gb = large_disk
    }
  }

  if (!is_small_callset) {
    # For large callsets we still need to gather the sharded metrics.
    call Tasks.GatherVariantCallingMetrics {
      input:
        input_details = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = callset_name,
        disk_size_gb = medium_disk
    }
  }

  # CrossCheckFingerprints takes forever on large callsets.
  # We scatter over the input GVCFs to make things faster.
  if (scatter_cross_check_fingerprints) {
    call Tasks.GetFingerprintingIntervalIndices {
      input:
        unpadded_intervals = unpadded_intervals,
        haplotype_database = haplotype_database
    }

    Array[Int] fingerprinting_indices = GetFingerprintingIntervalIndices.indices_to_fingerprint

    scatter (idx in fingerprinting_indices) {
      File vcfs_to_fingerprint = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx]
    }

    call Tasks.GatherVcfs as GatherFingerprintingVcfs {
      input:
        input_vcfs = vcfs_to_fingerprint,
        output_vcf_name = callset_name + ".gathered.fingerprinting.vcf.gz",
        disk_size_gb = medium_disk
    }

    call Tasks.SelectFingerprintSiteVariants {
      input:
        input_vcf = GatherFingerprintingVcfs.output_vcf,
        base_output_name = callset_name + ".fingerprinting",
        haplotype_database = haplotype_database,
        disk_size_gb = medium_disk
    }

    call Tasks.PartitionSampleNameMap {
      input:
        sample_name_map = sample_name_map,
        line_limit = 1000
    }

    scatter (idx in range(length(PartitionSampleNameMap.partitions))) {

      Array[File] files_in_partition = read_lines(PartitionSampleNameMap.partitions[idx])

      call Tasks.CrossCheckFingerprint as CrossCheckFingerprintsScattered {
        input:
          gvcf_paths = files_in_partition,
          vcf_paths = vcfs_to_fingerprint,
          sample_name_map = sample_name_map,
          haplotype_database = haplotype_database,
          output_base_name = callset_name + "." + idx,
          scattered = true
      }
    }

    call Tasks.GatherPicardMetrics as GatherFingerprintingMetrics {
      input:
        metrics_files = CrossCheckFingerprintsScattered.crosscheck_metrics,
        output_file_name = callset_name + ".fingerprintcheck",
        disk_size_gb = small_disk
    }
  }

  if (!scatter_cross_check_fingerprints) {

    scatter (line in sample_name_map_lines) {
      File gvcf_paths = line[1]
    }

    call Tasks.CrossCheckFingerprint as CrossCheckFingerprintSolo {
      input:
        gvcf_paths = gvcf_paths,
        vcf_paths = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
        sample_name_map = sample_name_map,
        haplotype_database = haplotype_database,
        output_base_name = callset_name
    }
  }

  # Get the metrics from either code path
  File output_detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherVariantCallingMetrics.detail_metrics_file])
  File output_summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherVariantCallingMetrics.summary_metrics_file])

  # Get the VCFs from either code path
  Array[File?] output_vcf_files = if defined(FinalGatherVcf.output_vcf) then [FinalGatherVcf.output_vcf] else FindFilteringThresholdAndFilter.output_vcf
  Array[File?] output_vcf_index_files = if defined(FinalGatherVcf.output_vcf_index) then [FinalGatherVcf.output_vcf_index] else FindFilteringThresholdAndFilter.output_vcf_index

  output {
    # Metrics from either the small or large callset
    File detail_metrics_file = output_detail_metrics_file
    File summary_metrics_file = output_summary_metrics_file

    Array[File] output_vcfs = select_all(output_vcf_files)
    Array[File] output_vcf_indices = select_all(output_vcf_index_files)

    # Output the interval list generated/used by this run workflow.
    Array[File] output_intervals = SplitIntervalList.output_intervals

    # Output the metrics from crosschecking fingerprints.
    File crosscheck_fingerprint_check = select_first([CrossCheckFingerprintSolo.crosscheck_metrics, GatherFingerprintingMetrics.gathered_metrics])
  }
  meta {
    allowNestedInputs: true
  }
}