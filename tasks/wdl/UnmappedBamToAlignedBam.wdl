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

import "./Alignment.wdl" as Alignment
import "./DragmapAlignment.wdl" as DragmapAlignment
import "./SplitLargeReadGroup.wdl" as SplitRG
import "./Qc.wdl" as QC
import "./BamProcessing.wdl" as Processing
import "./Utilities.wdl" as Utils
import "../../structs/dna_seq/DNASeqStructs.wdl" as Structs

# WORKFLOW DEFINITION
workflow UnmappedBamToAlignedBam {

  input {
    SampleAndUnmappedBams sample_and_unmapped_bams
    DNASeqSingleSampleReferences references
    DragmapReference? dragmap_reference
    PapiSettings papi_settings

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String cross_check_fingerprints_by
    File haplotype_database_file
    Float lod_threshold
    String recalibrated_bam_basename
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
    Boolean bin_base_qualities = true
    Boolean somatic = false
    Boolean perform_bqsr = true
    Boolean use_bwa_mem = true
    Boolean allow_empty_ref_alt = false
  }

  Float cutoff_for_large_rg_in_gb = 20.0

  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  Int compression_level = 2

  # Get the size of the standard reference files as well as the additional reference files needed for BWA

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in sample_and_unmapped_bams.flowcell_unmapped_bams) {

    Float unmapped_bam_size = size(unmapped_bam, "GiB")

    String unmapped_bam_basename = basename(unmapped_bam, sample_and_unmapped_bams.unmapped_bam_suffix)

    # QC the unmapped BAM
    call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics {
      input:
        input_bam = unmapped_bam,
        metrics_filename = unmapped_bam_basename + ".unmapped.quality_yield_metrics",
        preemptible_tries = papi_settings.preemptible_tries
    }

    if (unmapped_bam_size > cutoff_for_large_rg_in_gb) {
      # Split bam into multiple smaller bams,
      # map reads to reference and recombine into one bam
      call SplitRG.SplitLargeReadGroup as SplitRG {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
          reference_fasta = references.reference_fasta,
          dragmap_reference = dragmap_reference,
          compression_level = compression_level,
          preemptible_tries = papi_settings.preemptible_tries,
          hard_clip_reads = hard_clip_reads,
          unmap_contaminant_reads = unmap_contaminant_reads,
          use_bwa_mem = use_bwa_mem,
          allow_empty_ref_alt = allow_empty_ref_alt
      }
    }

    if (unmapped_bam_size <= cutoff_for_large_rg_in_gb) {
      # Map reads to reference
      if (use_bwa_mem) {
        call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
          input:
            input_bam = unmapped_bam,
            bwa_commandline = bwa_commandline,
            output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
            reference_fasta = references.reference_fasta,
            compression_level = compression_level,
            preemptible_tries = papi_settings.preemptible_tries,
            hard_clip_reads = hard_clip_reads,
            unmap_contaminant_reads = unmap_contaminant_reads,
            allow_empty_ref_alt = allow_empty_ref_alt
        }
      }
      if (!use_bwa_mem) {
        call DragmapAlignment.SamToFastqAndDragmapAndMba as SamToFastqAndDragmapAndMba {
          input:
            input_bam = unmapped_bam,
            output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
            reference_fasta = references.reference_fasta,
            dragmap_reference = select_first([dragmap_reference]),
            compression_level = compression_level,
            preemptible_tries = papi_settings.preemptible_tries,
            hard_clip_reads = hard_clip_reads,
            unmap_contaminant_reads = unmap_contaminant_reads
        }
      }
    }

    File output_aligned_bam = select_first([SamToFastqAndBwaMemAndMba.output_bam, SamToFastqAndDragmapAndMba.output_bam, SplitRG.aligned_bam])

    Float mapped_bam_size = size(output_aligned_bam, "GiB")

    # QC the aligned but unsorted readgroup BAM
    # no reference as the input here is unsorted, providing a reference would cause an error
    call QC.CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics {
      input:
        input_bam = output_aligned_bam,
        output_bam_prefix = unmapped_bam_basename + ".readgroup",
        preemptible_tries = papi_settings.preemptible_tries
    }
  }

  # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
  Float gb_size_cutoff_for_preemptibles = 110.0
  Boolean data_too_large_for_preemptibles = size(output_aligned_bam, "GiB") > gb_size_cutoff_for_preemptibles

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call Processing.MarkDuplicates as MarkDuplicates {
    input:
      input_bams = output_aligned_bam,
      output_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_and_unmapped_bams.base_file_name + ".duplicate_metrics",
      total_input_size = size(output_aligned_bam, "GiB"),
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
  }

  # Sort aggregated+deduped BAM file and fix tags
  call Processing.SortSam as SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicate_marked.sorted",
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
  }

  Float agg_bam_size = size(SortSampleBam.output_bam, "GiB")

  if (defined(haplotype_database_file)) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [ SortSampleBam.output_bam ],
        input_bam_indexes = [SortSampleBam.output_bam_index],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = sample_and_unmapped_bams.base_file_name + ".crosscheck",
        total_input_size = agg_bam_size,
        lod_threshold = lod_threshold,
        cross_check_by = cross_check_fingerprints_by,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Create list of sequences for scatter-gather parallelization
  call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
      ref_dict = references.reference_fasta.ref_dict,
      preemptible_tries = papi_settings.preemptible_tries
  }

  # Estimate level of cross-sample contamination
  call Processing.CheckContamination as CheckContamination {
    input:
      input_bam = SortSampleBam.output_bam,
      input_bam_index = SortSampleBam.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      output_prefix = sample_and_unmapped_bams.base_file_name + ".preBqsr",
      preemptible_tries = papi_settings.agg_preemptible_tries,
      contamination_underestimation_factor = 0.75
  }

  # We need disk to localize the sharded input and output due to the scatter for BQSR.
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space
  # to account for the fact that the data is not split evenly.
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel

  if (perform_bqsr) {
    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
      # Generate the recalibration model by interval
      call Processing.BaseRecalibrator as BaseRecalibrator {
        input:
          input_bam = SortSampleBam.output_bam,
          input_bam_index = SortSampleBam.output_bam_index,
          recalibration_report_filename = sample_and_unmapped_bams.base_file_name + ".recal_data.csv",
          sequence_group_interval = subgroup,
          dbsnp_vcf = references.dbsnp_vcf,
          dbsnp_vcf_index = references.dbsnp_vcf_index,
          known_indels_sites_vcfs = references.known_indels_sites_vcfs,
          known_indels_sites_indices = references.known_indels_sites_indices,
          ref_dict = references.reference_fasta.ref_dict,
          ref_fasta = references.reference_fasta.ref_fasta,
          ref_fasta_index = references.reference_fasta.ref_fasta_index,
          bqsr_scatter = bqsr_divisor,
          preemptible_tries = papi_settings.agg_preemptible_tries
      }
    }

    # Merge the recalibration reports resulting from by-interval recalibration
    # The reports are always the same size
    call Processing.GatherBqsrReports as GatherBqsrReports {
      input:
        input_bqsr_reports = BaseRecalibrator.recalibration_report,
        output_report_filename = sample_and_unmapped_bams.base_file_name + ".recal_data.csv",
        preemptible_tries = papi_settings.preemptible_tries
    }

    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
      # Apply the recalibration model by interval
      call Processing.ApplyBQSR as ApplyBQSR {
        input:
          input_bam = SortSampleBam.output_bam,
          input_bam_index = SortSampleBam.output_bam_index,
          output_bam_basename = recalibrated_bam_basename,
          recalibration_report = GatherBqsrReports.output_bqsr_report,
          sequence_group_interval = subgroup,
          ref_dict = references.reference_fasta.ref_dict,
          ref_fasta = references.reference_fasta.ref_fasta,
          ref_fasta_index = references.reference_fasta.ref_fasta_index,
          bqsr_scatter = bqsr_divisor,
          compression_level = compression_level,
          preemptible_tries = papi_settings.agg_preemptible_tries,
          bin_base_qualities = bin_base_qualities,
          somatic = somatic
      }
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call Processing.GatherSortedBamFiles as GatherBamFiles {
    input:
      input_bams = select_first([ApplyBQSR.recalibrated_bam, [SortSampleBam.output_bam]]),
      output_bam_basename = sample_and_unmapped_bams.base_file_name,
      total_input_size = agg_bam_size,
      compression_level = compression_level,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics

    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File? output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File output_bam = GatherBamFiles.output_bam
    File output_bam_index = GatherBamFiles.output_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}
