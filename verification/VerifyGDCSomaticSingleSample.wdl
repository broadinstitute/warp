version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as Tasks

workflow VerifyGDCSomaticSingleSample {

  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_bam
    File truth_bai
    File test_bam
    File test_bai
  }

  call MetricsVerification.VerifyMetrics as CompareMetrics {
    input:
      test_metrics = test_metrics,
      truth_metrics = truth_metrics
  }

  call Tasks.CompareBams {
    input:
      test_bam = test_bam,
      truth_bam = truth_bam,
  }

  output {
    Array[File] metric_comparison_report_files = CompareMetrics.metric_comparison_report_files
  }
  meta {
    allowNestedInputes: true
  }
}

