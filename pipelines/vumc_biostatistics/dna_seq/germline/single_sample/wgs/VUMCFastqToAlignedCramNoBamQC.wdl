version 1.0

## Copyright Broad Institute/VUMC, 2018/2023
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
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"
import "VUMCUnmappedBamToAlignedCramNoBamQC.wdl" as ToAlignedCram

# WORKFLOW DEFINITION
workflow VUMCFastqToAlignedCramNoBamQC {


  String pipeline_version = "3.1.7.1.beta"


  input {
    # Optional for VUMC pipeline
    String sample_name 
    String fastq_1 
    String fastq_2 
    String readgroup_name 
    String? library_name 
    String? platform_unit 
    String? run_date 
    String? platform_name 
    String? sequencing_center 

    # Optional for BROAD pipeline
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
    flowcell_unmapped_bams: [ PairedFastQsToUnmappedBAM.output_unmapped_bam ],
    sample_name: sample_name,
    unmapped_bam_suffix: ".bam"
  }

  call ToAlignedCram.VUMCUnmappedBamToAlignedCramNoBamQC as ToCram {
    input:
      sample_and_unmapped_bams = sample_and_unmapped_bams,
      references = references,
      dragmap_reference = dragmap_reference,
      papi_settings = papi_settings,

      dragen_functional_equivalence_mode = dragen_functional_equivalence_mode,
      dragen_maximum_quality_mode = dragen_maximum_quality_mode,

      unmap_contaminant_reads = unmap_contaminant_reads,
      perform_bqsr = perform_bqsr,
      use_bwa_mem = use_bwa_mem,
      allow_empty_ref_alt = allow_empty_ref_alt
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = ToCram.quality_yield_metrics

    File duplicate_metrics = ToCram.duplicate_metrics
    File? output_bqsr_reports = ToCram.output_bqsr_reports

    File output_cram = ToCram.output_cram
    File output_cram_index = ToCram.output_cram_index
    File output_cram_md5 = ToCram.output_cram_md5
  }
  meta {
    allowNestedInputs: true
  }
}
