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

import "../../../../../../tasks/vumc_biostatistics/VUMCUnmappedBamToAlignedBamNoBamQC.wdl" as ToBam
import "../../../../../../tasks/broad/AggregatedBamQC.wdl" as AggregatedQC
import "../../../../../../tasks/broad/Utilities.wdl" as Utilities

## Important notes by Quanhu Sheng, 20230815
## Higly recommended to use the following preset arguments:
## ConvertToCram:
##   disk_size = 100
## GatherBamFiles:
##   additional_disk = 100
## MarkDuplicates:
##   additional_disk = 200
##   memory_multiplier = 2
##   sorting_collection_size_ratio = 0.05
## SortSampleBam:
##   additional_disk = 100

workflow VUMCUnmappedBamToAlignedCramNoBamQC {
  input {
    SampleAndUnmappedBams sample_and_unmapped_bams
    DNASeqSingleSampleReferences references
    DragmapReference? dragmap_reference
    PapiSettings papi_settings

    Boolean dragen_functional_equivalence_mode = false
    Boolean dragen_maximum_quality_mode = false

    Boolean unmap_contaminant_reads = true
    Boolean perform_bqsr = true
    Boolean use_bwa_mem = true
    Boolean allow_empty_ref_alt = true
  }

  if (dragen_functional_equivalence_mode && dragen_maximum_quality_mode) {
    call Utilities.ErrorWithMessage as PresetArgumentsError {
      input:
        message = "Both dragen_functional_equivalence_mode and dragen_maximum_quality_mode have been set to true, however, they are mutually exclusive. You can set either of them to true, or set them both to false and adjust the arguments individually."
    }
  }

  # Set DRAGEN-related arguments according to the preset arguments
  Boolean unmap_contaminant_reads_ = if dragen_functional_equivalence_mode then false else (if dragen_maximum_quality_mode then true else unmap_contaminant_reads)
  Boolean perform_bqsr_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else perform_bqsr
  Boolean use_bwa_mem_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else use_bwa_mem

  # Not overridable:
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicates_marked.recalibrated"

  call ToBam.VUMCUnmappedBamToAlignedBamNoBamQC as UnmappedBamToAlignedBam {
    input:
      sample_and_unmapped_bams    = sample_and_unmapped_bams,
      references                  = references,
      dragmap_reference           = dragmap_reference,
      papi_settings               = papi_settings,

      lod_threshold               = lod_threshold,
      recalibrated_bam_basename   = recalibrated_bam_basename,
      perform_bqsr                = perform_bqsr_,
      use_bwa_mem                 = use_bwa_mem_,
      unmap_contaminant_reads     = unmap_contaminant_reads_,
      allow_empty_ref_alt         = allow_empty_ref_alt
  }

  call Utilities.ConvertToCram as ConvertToCram {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      output_basename = sample_and_unmapped_bams.base_file_name,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = UnmappedBamToAlignedBam.quality_yield_metrics

    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    File? output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
    File output_cram_md5 = ConvertToCram.output_cram_md5
  }
  meta {
    allowNestedInputs: true
  }
}
