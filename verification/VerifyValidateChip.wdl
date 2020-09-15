version 1.0

import "../verification/VerifyTasks.wdl" as Tasks

## Copyright Broad Institute, 2018
##
## This WDL script is designed to verify (compare) the outputs of a ValidateChip wdl.
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

workflow VerifyValidateChip {
  input {
    File test_gtc
    File truth_gtc
    File bead_pool_manifest_file

    File test_vcf
    File truth_vcf

    File test_genotype_concordance_vcf
    File truth_genotype_concordance_vcf

    File test_indel_genotype_concordance_vcf
    File truth_indel_genotype_concordance_vcf

    Array[File] test_metrics
    Array[File] truth_metrics
  }

  call Tasks.CompareGtcs {
    input:
      file1=test_gtc,
      file2=truth_gtc,
      bead_pool_manifest_file=bead_pool_manifest_file
  }

  call Tasks.CompareVcfs as CompareOutputVcfs {
    input:
      file1=truth_vcf,
      file2=test_vcf
  }

  call Tasks.CompareVcfs as CompareGenotypeConcordanceVcfs {
    input:
      file1=truth_genotype_concordance_vcf,
      file2=test_genotype_concordance_vcf
  }

  call Tasks.CompareVcfs as CompareIndelGenotypeConcordanceVcfs {
    input:
      file1=truth_indel_genotype_concordance_vcf,
      file2=test_indel_genotype_concordance_vcf
  }
  meta {
    allowNestedInputs: true
  }
}