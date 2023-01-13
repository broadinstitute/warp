version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as Tasks

## Copyright Broad Institute, 2021
##
## This WDL script is designed to verify (compare) the outputs of the CheckFingerprint wdl.
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

workflow VerifyCheckFingerprint {

  input {
    Array[File]? test_metrics
    Array[File]? truth_metrics

    File test_fingerprint_vcf
    File truth_fingerprint_vcf

    Boolean? done
  }

  if (defined(test_metrics) && defined(truth_metrics)) {
    call MetricsVerification.VerifyMetrics as CompareMetrics {
      input:
        test_metrics = select_first([test_metrics]),
        truth_metrics = select_first([truth_metrics])
    }
  }

  call Tasks.CompareVcfsAllowingQualityDifferences as CompareOutputFingerprintVcfs {
    input:
      file1 = test_fingerprint_vcf,
      file2 = truth_fingerprint_vcf
  }

  output {
    Array[File]? metric_comparison_report_files = CompareMetrics.metric_comparison_report_files
  }
  meta {
    allowNestedInputs: true
  }
}

