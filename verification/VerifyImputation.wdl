version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as Tasks

## Copyright Broad Institute, 2018
##
## This WDL script is designed to verify (compare) the outputs of an ArrayWf wdl.
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

workflow VerifyImputation {
  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_vcf
    File test_vcf

    Array[File]? single_sample_truth_vcf
    Array[File]? single_sample_test_vcf
  }

  # commenting out for now because failing on header,
  # might  need  a new program to compare these
  # metrics are not picard metrics
  # might consider diff because they are small files
  # call MetricsVerification.VerifyMetrics as CompareMetrics {
  #   input:
  #     test_metrics = test_metrics,
  #     truth_metrics = truth_metrics
  # }

  call Tasks.CompareVcfs as CompareOutputVcfs {
    input:
      file1 = truth_vcf,
      file2 = test_vcf
  }

  if (defined(single_sample_truth_vcf)) {
    scatter (idx in range(length(select_first([single_sample_truth_vcf])))) {
      call Tasks.CompareVcfs as CompareSingleSampleOutputVcfs {
        input:
          file1 = select_first([single_sample_test_vcf])[idx],
          file2 = select_first([single_sample_truth_vcf])[idx]
      }
    }
  }

  # TODO add fingerprint checks for verification

  output {
  }
  meta {
    allowNestedInputs: true
  }
}
