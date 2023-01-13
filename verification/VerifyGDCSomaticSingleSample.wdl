version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyGDCSomaticSingleSample {

  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_bam
    File truth_bai
    File test_bam
    File test_bai

    Boolean? done
  }

  call MetricsVerification.VerifyMetrics as CompareMetrics {
    input:
      test_metrics = test_metrics,
      truth_metrics = truth_metrics
  }

  call VerifyTasks.CompareBams {
    input:
      test_bam = test_bam,
      truth_bam = truth_bam,
      lenient_header = true
  }

  output {
    Array[File] metric_comparison_report_files = CompareMetrics.metric_comparison_report_files
  }
  meta {
    allowNestedInputes: true
  }
}

