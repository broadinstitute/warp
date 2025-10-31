version 1.0

# JointGenotypingByChromosomePartOne is now deprecated 2025-03-06

import "../../../../../../tasks/wdl/JointGenotypingTasks.wdl" as Tasks

# Joint Genotyping for hg38 Exomes and Whole Genomes (has not been tested on hg19)
workflow JointGenotypingByChromosomePartOne {

  String pipeline_version = "1.5.3"

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

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69

    File haplotype_database

    Int? top_level_scatter_count
    Float unbounded_scatter_count_scale_factor = 0.15
    Int gnarly_scatter_count = 10
    Boolean use_gnarly_genotyper = false

    # These workarounds exist to work around bugs in GenomicsDB.
    # By manually computing different shards, we can provide them using
    # The substitution mechanism provided by these workarounds.
    # The workaround works by providing the unique "callset_name.interval_name"
    # that failed as the key, and putting the manually-generated data as the
    # value in the `workaround` map. The key also needs to be added to the
    # workaround_keys` string, because WDL doesn't have a `Map.contains(...)`.
    # Keys should be concatenated together with '|' used as a separator.
    Map[String, File] gnarly_workaround = {"foo": sample_name_map}
    Map[String, File] gnarly_workaround_annotations = {"foo": sample_name_map}
    String gnarly_workaround_keys = "|"
    Map[String, File] import_workaround = {"foo": sample_name_map}
    String import_workaround_keys = "|"
  }

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000

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

    String import_workaround_key = callset_name + "." + basename(unpadded_intervals[idx])
    Boolean import_workaround_exists = sub(import_workaround_keys, import_workaround_key, "") != import_workaround_keys

    if (import_workaround_exists) {
      File faked_genomicsdb = import_workaround[import_workaround_key]
    }

    if (!import_workaround_exists) {
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
    }

    File genomicsdb = select_first([faked_genomicsdb, ImportGVCFs.output_genomicsdb])

    if (use_gnarly_genotyper) {

      # Further split the intervals for GnarlyGenotyper
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

        String gnarly_workaround_key = callset_name + "." + idx + "." + gnarly_idx
        Boolean gnarly_workaround_exists = sub(gnarly_workaround_keys, gnarly_workaround_key, "") != gnarly_workaround_keys

        if (gnarly_workaround_exists) {
          File faked_gnarly_genotyped_vcf = gnarly_workaround[gnarly_workaround_key]
        }

        if (!gnarly_workaround_exists) {
          call Tasks.GnarlyGenotyper {
            input:
              workspace_tar = genomicsdb,
              interval = gnarly_intervals[gnarly_idx],
              output_vcf_filename = gnarly_workaround_key + ".vcf.gz",
              ref_fasta = ref_fasta,
              ref_fasta_index = ref_fasta_index,
              ref_dict = ref_dict,
              dbsnp_vcf = dbsnp_vcf,
          }
        }

        File gnarly_vcf = select_first([faked_gnarly_genotyped_vcf, GnarlyGenotyper.output_vcf])
      }

      Array[File] gnarly_vcfs = gnarly_vcf

      # Gather the VCFs generated in the GnarlyGenotyper shards so that we end
      # up with the same number of VCF shards as the original top-level scatter
      call Tasks.GatherVcfs as TotallyRadicalGatherVcfs {
        input:
          input_vcfs = gnarly_vcfs,
          output_vcf_name = callset_name + "." + idx + ".gnarly.vcf.gz",
          disk_size_gb = large_disk
      }
  }


    if (!use_gnarly_genotyper) {
      call Tasks.GenotypeGVCFs {
        input:
          workspace_tar = genomicsdb,
          interval = unpadded_intervals[idx],
          output_vcf_filename = callset_name + "." + idx + "annotations.vcf.gz",
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
        disk_size_gb = medium_disk
    }
  }

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size_gb = medium_disk
  }


  call Tasks.GetFingerprintingIntervalIndices {
    input:
      unpadded_intervals = unpadded_intervals,
      haplotype_database = haplotype_database
  }

  # The result of GetFingerprintingIntervalIndices with no indices is []
  if (length(GetFingerprintingIntervalIndices.indices_to_fingerprint) > 0) {
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
  }

  output {
    # Outputs from the small callset path through the wdl.
    File output_vcf = SitesOnlyGatherVcf.output_vcf
    File output_vcf_index = SitesOnlyGatherVcf.output_vcf_index

    Array[File] output_hard_filtered_vcfs = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf
    Array[File] output_hard_filtered_vcf_indices = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index

    File? output_fingerprinting_vcf = SelectFingerprintSiteVariants.output_vcf
    File? output_fingerprinting_vcf_index = SelectFingerprintSiteVariants.output_vcf_index

    Array[File] genomics_databases = genomicsdb

    # Output the interval list generated/used by this run workflow.
    Array[File] output_intervals = SplitIntervalList.output_intervals
  }
  meta {
    allowNestedInputs: true
  }
}
