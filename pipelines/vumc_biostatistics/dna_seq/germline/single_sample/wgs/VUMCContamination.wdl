version 1.0

## Copyright Broad Institute/VUMC, 2018/2022
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

import "../../../../../../tasks/broad/BamProcessing.wdl" as Processing
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow VUMCContamination {

  String pipeline_version = "1.0.0.0"

  input {
    File input_cram
    File input_cram_index

    String sample_name

    DNASeqSingleSampleReferences references
    PapiSettings papi_settings
  }

  # Estimate level of cross-sample contamination
  call Processing.CheckContamination as CheckContamination {
    input:
      input_bam = input_cram,
      input_bam_index = input_cram_index,
      contamination_sites_ud = references.contamination_sites_ud,
      contamination_sites_bed = references.contamination_sites_bed,
      contamination_sites_mu = references.contamination_sites_mu,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      output_prefix = sample_name,
      preemptible_tries = papi_settings.agg_preemptible_tries,
      contamination_underestimation_factor = 0.75
  }

  # Outputs that will be retained when execution is complete
  output {
    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination
  }
  meta {
    allowNestedInputs: true
  }
}
