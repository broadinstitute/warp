version 1.0

import "../verification/VerifyTasks.wdl" as Tasks
import "../verification/VerifyMetrics.wdl" as MetricsVerification

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

    Boolean? done
  }

  call Tasks.CompareGtcs {
    input:
      file1=test_gtc,
      file2=truth_gtc,
      bead_pool_manifest_file=bead_pool_manifest_file
  }

  call Tasks.CompareVcfsAllowingQualityDifferences as CompareOutputVcfs {
    input:
      file1=truth_vcf,
      file2=test_vcf
  }

  call Tasks.CompareVcfsAllowingQualityDifferences as CompareGenotypeConcordanceVcfs {
    input:
      file1=truth_genotype_concordance_vcf,
      file2=test_genotype_concordance_vcf
  }

  call Tasks.CompareVcfsAllowingQualityDifferences as CompareIndelGenotypeConcordanceVcfs {
    input:
      file1=truth_indel_genotype_concordance_vcf,
      file2=test_indel_genotype_concordance_vcf
  }

  call MetricsVerification.CompareTwoNumbers {
    input:
      num1 = length(test_metrics),
      num2 = length(truth_metrics),
      error_msg = "Different number of metric files"
  }

  String avcdm_ext = "arrays_variant_calling_detail_metrics"

  scatter (idx in range(length(truth_metrics))) {
    String metrics_basename = basename(truth_metrics[idx])
    Boolean is_avcdm_file = basename(metrics_basename, avcdm_ext) != metrics_basename
    Array[String] metrics_to_ignore = if (is_avcdm_file) then ["AUTOCALL_DATE", "ANALYSIS_VERSION", "PIPELINE_VERSION"] else []
    call MetricsVerification.CompareMetricFiles {
      input:
        dependency_input = CompareTwoNumbers.output_file,
        file1 = test_metrics[idx],
        file2 = truth_metrics[idx],
        output_file = "metric_~{idx}.txt",
        metrics_to_ignore = metrics_to_ignore
    }
  }

  meta {
    allowNestedInputs: true
  }

  output {}
}