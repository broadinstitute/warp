version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as Tasks

workflow VerifyUltimaGenomicsWholeGenomeCramOnly {

  input {
    Array[File] test_metrics
    Array[File] truth_metrics

    File test_cram
    File truth_cram
    File test_crai
    File truth_crai

    Boolean? done
  }

  ## TEMP - Can't use this simple call, because the picard version used by UltimaGenomicsWholeGenomeGermline is based on an old version of Picard,
#  call MetricsVerification.VerifyMetrics as CompareMetrics {
#    input:
#      test_metrics = test_metrics,
#      truth_metrics = truth_metrics
#  }
  # So, I have to do all this - when the UltimaGenomicsWholeGenomeGermline picard is properly rebased, this should all be able to be rolled back.
  call MetricsVerification.CompareTwoNumbers {
    input:
      num1 = length(test_metrics),
      num2 = length(truth_metrics),
      error_msg = "Different number of metric files"
  }

  String alignment_metrics_ext = ".alignment_summary_metrics"
  String duplicate_metrics_ext = ".duplicate_metrics"

  scatter (idx in range(length(truth_metrics))) {
    String metrics_basename = basename(truth_metrics[idx])
    Boolean is_alignment_metrics_file = basename(metrics_basename, alignment_metrics_ext) != metrics_basename
    Boolean is_duplicate_metrics_file = basename(metrics_basename, duplicate_metrics_ext) != metrics_basename
    if ((!is_alignment_metrics_file) && (!is_duplicate_metrics_file)) {
      call MetricsVerification.CompareMetricFiles {
        input:
          dependency_input = CompareTwoNumbers.output_file,
          file1 = test_metrics[idx],
          file2 = truth_metrics[idx],
          output_file = "metric_~{idx}.txt",
          metrics_to_ignore = []
      }
    }
    if ((is_alignment_metrics_file) || (is_duplicate_metrics_file)) {
      call CompareOldMetricFiles {
        input:
          dependency_input = CompareTwoNumbers.output_file,
          file1 = test_metrics[idx],
          file2 = truth_metrics[idx],
          output_file = "metric_~{idx}.txt",
          metrics_to_ignore = []
      }
    }
  }
  ## END TEMP

  call Tasks.CompareCrams {
    input:
      test_cram = test_cram,
      test_crai = test_crai,
      truth_cram = truth_cram,
      truth_crai = truth_crai
  }

  call Tasks.CompareCrais {
    input:
      test_crai = test_crai,
      truth_crai = truth_crai
  }

  meta {
    allowNestedInputs: true
  }

  output{}
}



task CompareOldMetricFiles {
  input {
    File? dependency_input
    File file1
    File file2
    String output_file
    Array[String] metrics_to_ignore
  }

  command <<<
    java -Xmx3g -Dpicard.useLegacyParser=false  -jar /usr/gitc/picard.jar \
    CompareMetrics \
    --INPUT ~{file1} \
    --INPUT ~{file2} \
    --OUTPUT ~{output_file} \
    ~{true="--METRICS_TO_IGNORE" false="" length(metrics_to_ignore) > 0} ~{default="" sep=" --METRICS_TO_IGNORE " metrics_to_ignore}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:samtools_picard_bwa_snapshot_UG"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
  }
  output {
    File report_file = output_file
  }
}
