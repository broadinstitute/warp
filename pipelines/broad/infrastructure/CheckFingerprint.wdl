version 1.0

import "../../../tasks/broad/Utilities.wdl" as utils
import "../../../tasks/broad/InternalTasks.wdl" as InternalTasks
import "../../../tasks/broad/IlluminaGenotypingArrayTasks.wdl" as GenotypingTasks


## Copyright Broad Institute, 2021
##
## This WDL pipeline implements A CheckFingerprint Task
##
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

workflow CheckFingerprint {

  String pipeline_version = "0.1.0"

  input {
    # The name of the sample in the input_vcf.  Not required if there is only one sample in the VCF
    String? input_sample_alias
    File? input_vcf
    File? input_vcf_index
    File? input_bam
    File? input_bam_index

    String sample_alias
    String sample_lsid

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # If this is true, we will read fingerprints from Mercury
    # Otherwise, we will use the optional input fingerprint VCFs below
    Boolean read_fingerprint_from_mercury = false
    File? fingerprint_genotypes_vcf
    File? fingerprint_genotypes_vcf_index
    File haplotype_database_file

    String environment
    File vault_token_path
  }

  if (!defined(input_vcf) && !defined(input_bam)) {
    call utils.ErrorWithMessage as ErrorMessageInput {
      input:
        message = "Either input_vcf or input_bam (and NOT both) must be defined as input"
    }
  }

  if (defined(input_vcf) && defined(input_bam)) {
    call utils.ErrorWithMessage as ErrorMessageDoubleInput {
      input:
        message = "input_vcf and input_bam cannot both be defined as input"
    }
  }

  call InternalTasks.MakeSafeFilename {
    input:
      name = sample_alias
  }

  if (read_fingerprint_from_mercury) {
    call InternalTasks.DownloadGenotypes {
      input:
        sample_alias = sample_alias,
        sample_lsid = sample_lsid,
        output_vcf_base_name = MakeSafeFilename.output_safe_name + ".reference.fingerprint",
        haplotype_database_file = haplotype_database_file,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        environment = environment,
        vault_token_path = vault_token_path
    }
  }

  Boolean fingerprint_downloaded_from_mercury = select_first([DownloadGenotypes.fingerprint_retrieved, false])

  File? fingerprint_vcf_to_use = if (fingerprint_downloaded_from_mercury) then DownloadGenotypes.reference_fingerprint_vcf else fingerprint_genotypes_vcf
  File? fingerprint_vcf_index_to_use = if (fingerprint_downloaded_from_mercury) then DownloadGenotypes.reference_fingerprint_vcf_index else fingerprint_genotypes_vcf_index

  if (defined(fingerprint_vcf_to_use)) {
    call GenotypingTasks.CheckFingerprint {
      input:
        input_vcf_file = input_vcf,
        input_vcf_index_file = input_vcf_index,
        input_sample_alias = input_sample_alias,
        genotypes_vcf_file = select_first([fingerprint_vcf_to_use]),
        genotypes_vcf_index_file = select_first([fingerprint_vcf_to_use]),
        haplotype_database_file = haplotype_database_file,
        expected_sample_alias = sample_alias,
        output_metrics_basename = MakeSafeFilename.output_safe_name
    }
  }

  output {
    File? reference_fingerprint_vcf = DownloadGenotypes.reference_fingerprint_vcf
    File? reference_fingerprint_vcf_index = DownloadGenotypes.reference_fingerprint_vcf_index
    File? fingerprint_detail_metrics_file = CheckFingerprint.detail_metrics
    File? fingerprint_summary_metrics_file = CheckFingerprint.summary_metrics
    Float? lod_score = CheckFingerprint.lod
  }
  meta {
    allowNestedInputs: true
  }
}
