version 1.0

import "../../../../../tasks/wdl/JointGenotypingTasks.wdl" as Tasks
import "https://raw.githubusercontent.com/broadinstitute/gatk/4.5.0.0/scripts/vcf_site_level_filtering_wdl/JointVcfFiltering.wdl" as Filtering


# Joint Genotyping for hg38 Whole Genomes and Exomes (has not been tested on hg19)
workflow JointGenotyping {

  String pipeline_version = "1.7.3"

  input {
    File unpadded_intervals_file

    String callset_name
    File sample_name_map

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    Int small_disk
    Int medium_disk
    Int large_disk
    Int huge_disk

    Array[String]? snp_recalibration_tranche_values
    Array[String] snp_recalibration_annotation_values
    Array[String]? indel_recalibration_tranche_values
    Array[String]? indel_recalibration_annotation_values

    File haplotype_database

    File eval_interval_list
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf = dbsnp_vcf
    File dbsnp_resource_vcf_index = dbsnp_vcf_index

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69
    Float? vqsr_snp_filter_level
    Float? vqsr_indel_filter_level
    Int? snp_vqsr_downsampleFactor
    File? targets_interval_list

    Int? top_level_scatter_count
    Boolean? gather_vcfs
    Int snps_variant_recalibration_threshold = 500000
    Boolean rename_gvcf_samples = true
    Float unbounded_scatter_count_scale_factor = 0.15
    Int gnarly_scatter_count = 10
    Boolean use_gnarly_genotyper = false
    Boolean use_allele_specific_annotations = true # only applicabale to VQSR
    Boolean cross_check_fingerprints = true
    Boolean scatter_cross_check_fingerprints = false
    Boolean run_vets = false
  }

  Boolean allele_specific_annotations = !use_gnarly_genotyper && use_allele_specific_annotations

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

    if (use_gnarly_genotyper) {

      call Tasks.SplitIntervalList as GnarlyIntervalScatterDude {
        input:
          interval_list = unpadded_intervals[idx],
          scatter_count = gnarly_scatter_count,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size_gb = small_disk,
          sample_names_unique_done = CheckSamplesUnique.samples_unique
      }

      Array[File] gnarly_intervals = GnarlyIntervalScatterDude.output_intervals

      scatter (gnarly_idx in range(length(gnarly_intervals))) {
        call Tasks.GnarlyGenotyper {
          input:
            workspace_tar = ImportGVCFs.output_genomicsdb,
            interval = gnarly_intervals[gnarly_idx],
            output_vcf_filename = callset_name + "." + idx + "." + gnarly_idx + ".vcf.gz",
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp_vcf = dbsnp_vcf,
        }
      }

      Array[File] gnarly_gvcfs = GnarlyGenotyper.output_vcf

      call Tasks.GatherVcfs as TotallyRadicalGatherVcfs {
        input:
          input_vcfs = gnarly_gvcfs,
          output_vcf_name = callset_name + "." + idx + ".gnarly.vcf.gz",
          disk_size_gb = large_disk
      }
    }

    if (!use_gnarly_genotyper) {
      call Tasks.GenotypeGVCFs {
        input:
          workspace_tar = ImportGVCFs.output_genomicsdb,
          interval = unpadded_intervals[idx],
          output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf,
          disk_size_gb = medium_disk
      }
    }

    File genotyped_vcf = select_first([TotallyRadicalGatherVcfs.output_vcf, GenotypeGVCFs.output_vcf])
    File genotyped_vcf_index = select_first([TotallyRadicalGatherVcfs.output_vcf_index, GenotypeGVCFs.output_vcf_index])

    call Tasks.HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = genotyped_vcf,
        vcf_index = genotyped_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        targets_interval_list = targets_interval_list,
        disk_size_gb = medium_disk
    }
  }

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size_gb = medium_disk
  }
  
  if (run_vets) {
    String resource_args = " --resource:hapmap,training=true,calibration=true " + hapmap_resource_vcf + 
      " --resource:omni,training=true,calibration=true " + omni_resource_vcf + 
      " --resource:1000G,training=true " + one_thousand_genomes_resource_vcf +
      " --resource:mills,training=true,calibration=true " + mills_resource_vcf + " "
    String extract_extra_args = if defined(targets_interval_list) then " -L " + targets_interval_list + " " else "" #only train the model over the targets, apply the model to everything

    call Filtering.JointVcfFiltering as TrainAndApplyVETS {
      input:
        input_vcfs = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
        input_vcf_idxs = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index,
        sites_only_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_vcf_idx = SitesOnlyGatherVcf.output_vcf_index,
        annotations = snp_recalibration_annotation_values, #the snp list here is a superset of  the indel list
        extract_extra_args = extract_extra_args,
        resource_args = resource_args,
        output_prefix = callset_name,
        gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
  } 
  if (!run_vets) {
    call Tasks.IndelsVariantRecalibrator {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".indels.recal",
        tranches_filename = callset_name + ".indels.tranches",
        recalibration_tranche_values = select_first([indel_recalibration_tranche_values]),
        recalibration_annotation_values = select_first([indel_recalibration_annotation_values]),
        mills_resource_vcf = mills_resource_vcf,
        mills_resource_vcf_index = mills_resource_vcf_index,
        axiomPoly_resource_vcf = axiomPoly_resource_vcf,
        axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size_gb = small_disk
    }

    if (num_gvcfs > snps_variant_recalibration_threshold) {
      call Tasks.SNPsVariantRecalibratorCreateModel {
        input:
          sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
          sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
          recalibration_filename = callset_name + ".snps.recal",
          tranches_filename = callset_name + ".snps.tranches",
          recalibration_tranche_values = select_first([snp_recalibration_tranche_values]),
          recalibration_annotation_values = snp_recalibration_annotation_values,
          downsampleFactor = select_first([snp_vqsr_downsampleFactor]),
          model_report_filename = callset_name + ".snps.model.report",
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
          use_allele_specific_annotations = allele_specific_annotations,
          disk_size_gb = small_disk
      }

      scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {
        call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
          input:
            sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
            sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
            recalibration_filename = callset_name + ".snps." + idx + ".recal",
            tranches_filename = callset_name + ".snps." + idx + ".tranches",
            recalibration_tranche_values = select_first([snp_recalibration_tranche_values]),
            recalibration_annotation_values = snp_recalibration_annotation_values,
            model_report = SNPsVariantRecalibratorCreateModel.model_report,
            hapmap_resource_vcf = hapmap_resource_vcf,
            hapmap_resource_vcf_index = hapmap_resource_vcf_index,
            omni_resource_vcf = omni_resource_vcf,
            omni_resource_vcf_index = omni_resource_vcf_index,
            one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
            one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
            dbsnp_resource_vcf = dbsnp_resource_vcf,
            dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
            use_allele_specific_annotations = allele_specific_annotations,
            disk_size_gb = small_disk
          }
      }

      call Tasks.GatherTranches as SNPGatherTranches {
        input:
          tranches = SNPsVariantRecalibratorScattered.tranches,
          output_filename = callset_name + ".snps.gathered.tranches",
          mode = "SNP",
          disk_size_gb = small_disk
      }
    }

    if (num_gvcfs <= snps_variant_recalibration_threshold) {
      call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
        input:
          sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
          sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
          recalibration_filename = callset_name + ".snps.recal",
          tranches_filename = callset_name + ".snps.tranches",
          recalibration_tranche_values = select_first([snp_recalibration_tranche_values]),
          recalibration_annotation_values = snp_recalibration_annotation_values,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
          use_allele_specific_annotations = allele_specific_annotations,
          disk_size_gb = small_disk
      }
    }

    scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
      #for really large callsets we give to friends, just apply filters to the sites-only
      call Tasks.ApplyRecalibration {
        input:
          recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
          input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
          input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
          indels_recalibration = IndelsVariantRecalibrator.recalibration,
          indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
          indels_tranches = IndelsVariantRecalibrator.tranches,
          snps_recalibration = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
          snps_recalibration_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
          snps_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
          indel_filter_level = select_first([vqsr_indel_filter_level]),
          snp_filter_level = select_first([vqsr_snp_filter_level]),
          use_allele_specific_annotations = allele_specific_annotations,
          disk_size_gb = medium_disk
      }
    }
  }

  Array[File] filtered_vcfs = select_first([ApplyRecalibration.recalibrated_vcf, TrainAndApplyVETS.scored_vcfs])
  Array[File] filtered_vcf_idxs = select_first([ApplyRecalibration.recalibrated_vcf_index, TrainAndApplyVETS.scored_vcf_idxs])

  if (!is_small_callset) {
    scatter (idx in range(length(filtered_vcfs))) {
    # For large callsets we need to collect metrics from the shards and gather them later.
      call Tasks.CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = filtered_vcfs[idx],
          input_vcf_index = filtered_vcf_idxs[idx],
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size_gb = medium_disk
      }
    }

    call Tasks.GatherVariantCallingMetrics {
      input:
        input_details = CollectMetricsSharded.detail_metrics_file,
        input_summaries = CollectMetricsSharded.summary_metrics_file,
        output_prefix = callset_name,
        disk_size_gb = medium_disk
    }
  }

  # For small callsets we can gather the VCF shards and then collect metrics on it.
  if (is_small_callset) {
    call Tasks.GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs = filtered_vcfs,
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

  # CrossCheckFingerprints takes forever on large callsets.
  # We scatter over the input GVCFs to make things faster.
if (cross_check_fingerprints) {
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
        vcf_paths = filtered_vcfs,
        sample_name_map = sample_name_map,
        haplotype_database = haplotype_database,
        output_base_name = callset_name
      }
    }
    File crosscheck_fingerprint_results = select_first([CrossCheckFingerprintSolo.crosscheck_metrics, GatherFingerprintingMetrics.gathered_metrics])
  }

  # Get the metrics from either code path
  File output_detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherVariantCallingMetrics.detail_metrics_file])
  File output_summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherVariantCallingMetrics.summary_metrics_file])

  # Get the VCFs from either code path
  Array[File?] output_vcf_files = if defined(FinalGatherVcf.output_vcf) then [FinalGatherVcf.output_vcf] else filtered_vcfs
  Array[File?] output_vcf_index_files = if defined(FinalGatherVcf.output_vcf_index) then [FinalGatherVcf.output_vcf_index] else filtered_vcf_idxs

  output {
    # Metrics from either the small or large callset
    File detail_metrics_file = output_detail_metrics_file
    File summary_metrics_file = output_summary_metrics_file

    # Outputs from the small callset path through the wdl.
    Array[File] output_vcfs = select_all(output_vcf_files)
    Array[File] output_vcf_indices = select_all(output_vcf_index_files)

    # Output the interval list generated/used by this run workflow.
    Array[File] output_intervals = SplitIntervalList.output_intervals

    # Output the metrics from crosschecking fingerprints.
    File? crosscheck_fingerprint_check = crosscheck_fingerprint_results
  }
  meta {
    allowNestedInputs: true
  }
}
