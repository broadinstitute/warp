version 1.0
# The JointGenotypingByChromosomePartTwo is deprecated 2025-03-06
import "../../../../../../tasks/wdl/JointGenotypingTasks.wdl" as Tasks

# Joint Genotyping for hg38 Exomes and Whole Genomes (has not been tested on hg19)
workflow JointGenotypingByChromosomePartTwo {

  String pipeline_version = "1.5.3"

  input {
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

    File sites_only_vcfs_fofn
    File sites_only_vcf_indices_fofn
    File hard_filtered_with_genotypes_vcfs_fofn
    File hard_filtered_with_genotypes_vcf_indices_fofn
    File? annotation_db_vcfs_fofn
    File? annotation_db_vcf_indices_fofn
    File fingerprinting_vcfs_fofn
    File fingerprinting_vcf_indices_fofn

    Array[String] snp_recalibration_tranche_values
    Array[String] snp_recalibration_annotation_values
    Array[String] indel_recalibration_tranche_values
    Array[String] indel_recalibration_annotation_values

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
    Float snp_filter_level
    Float indel_filter_level
    Int SNP_VQSR_downsampleFactor

    Boolean use_gnarly_genotyper = false
    Boolean use_allele_specific_annotations = true
  }

  Boolean allele_specific_annotations = !use_gnarly_genotyper && use_allele_specific_annotations
  Boolean make_annotation_db = false

  Array[File] sites_only_vcfs = read_lines(sites_only_vcfs_fofn)
  Array[File] sites_only_vcf_indices = read_lines(sites_only_vcf_indices_fofn)
  Array[File] hard_filtered_with_genotypes_vcfs = read_lines(hard_filtered_with_genotypes_vcfs_fofn)
  Array[File] hard_filtered_with_genotypes_vcf_indices = read_lines(hard_filtered_with_genotypes_vcf_indices_fofn)
  Array[File] fingerprinting_vcfs = read_lines(fingerprinting_vcfs_fofn)
  Array[File] fingerprinting_vcf_indices = read_lines(fingerprinting_vcf_indices_fofn)

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs = sites_only_vcfs,
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size_gb = medium_disk
  }

  if (use_gnarly_genotyper && make_annotation_db) {
    Array[File] annotation_db_vcfs = read_lines(select_first([annotation_db_vcfs_fofn, ""]))
    Array[File] annotation_db_vcf_indices = read_lines(select_first([annotation_db_vcf_indices_fofn, ""]))
    call Tasks.GatherVcfs as GatherAnnotationDBVcf {
      input:
        input_vcfs = annotation_db_vcfs,
        output_vcf_name = callset_name + ".annotationDB.vcf.gz",
        disk_size_gb = medium_disk
    }
  }

  call Tasks.IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = allele_specific_annotations,
      disk_size_gb = small_disk
  }

  call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".snps.recal",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_tranche_values = snp_recalibration_tranche_values,
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

  scatter (idx in range(length(sites_only_vcfs))) {
    #for super large callsets only apply recalibration to the sites-only
    call Tasks.ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = sites_only_vcfs[idx],
        input_vcf_index = sites_only_vcf_indices[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = SNPsVariantRecalibratorClassic.recalibration,
        snps_recalibration_index = SNPsVariantRecalibratorClassic.recalibration_index,
        snps_tranches = SNPsVariantRecalibratorClassic.tranches,
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size_gb = medium_disk
    }
  }

  scatter (idx in range(length(hard_filtered_with_genotypes_vcfs))) {
    # For large callsets we need to collect metrics from the shards and gather them later.
    call Tasks.CollectVariantCallingMetrics as CollectMetricsSharded {
      input:
        input_vcf = hard_filtered_with_genotypes_vcfs[idx],
        input_vcf_index = hard_filtered_with_genotypes_vcf_indices[idx],
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
        vcf_paths = fingerprinting_vcfs,
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

  output {
    File? annotation_db_vcf = GatherAnnotationDBVcf.output_vcf
    File? annotation_db_vcf_index = GatherAnnotationDBVcf.output_vcf_index

    # Metrics from either the small or large callset
    File detail_metrics_file = GatherVariantCallingMetrics.detail_metrics_file
    File summary_metrics_file = GatherVariantCallingMetrics.summary_metrics_file

    # Outputs from the small callset path through the wdl.
    Array[File] output_vcfs = ApplyRecalibration.recalibrated_vcf
    Array[File] output_vcf_indices = ApplyRecalibration.recalibrated_vcf_index

    File crosscheck_metrics = GatherFingerprintingMetrics.gathered_metrics
  }
  meta {
    allowNestedInputs: true
  }
}

