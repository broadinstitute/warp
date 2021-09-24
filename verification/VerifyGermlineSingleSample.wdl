version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as Tasks

workflow VerifyGermlineSingleSample {

  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_cram
    File truth_crai
    File test_cram
    File test_crai

    File truth_gvcf
    File test_gvcf
  }

  call MetricsVerification.VerifyMetrics as CompareMetrics {
    input:
      test_metrics = test_metrics,
      truth_metrics = truth_metrics
  }

  call CompareGvcfs {
    input:
      test_gvcf = test_gvcf,
      truth_gvcf = truth_gvcf
  }

  call Tasks.CompareCrais {
    input:
      test_crai = test_crai,
      truth_crai = truth_crai
  }

  call Tasks.CompareCrams {
    input:
      test_cram = test_cram,
      test_crai = test_crai,
      truth_cram = truth_cram,
      truth_crai = truth_crai
  }

  output {
    Array[File] metric_comparison_report_files = CompareMetrics.metric_comparison_report_files
  }
  meta {
    allowNestedInputs: true
  }
}

task CompareGvcfs {

  input {
    File test_gvcf
    File truth_gvcf

  }

  command {
    set -e -o pipefail

    exit_code=0

    if cmp <( gunzip -c -f ~{test_gvcf} | grep -v '^##' ) <( gunzip -c -f ~{truth_gvcf} | grep -v '^##' ); then
      exit 0
    fi

    /usr/bin/diff <( gunzip -c -f ~{test_gvcf} | grep -v '^##' ) <( gunzip -c -f ~{truth_gvcf} | grep -v '^##' ) > gvcf_diff.txt

    DIFF_LINES=$( grep -e "^<" gvcf_diff.txt | wc -l )

    echo "$DIFF_LINES" > diff_lines_total.txt

    if [ $DIFF_LINES -ge 10 ]; then
      exit_code=1
      
      echo "Error: GVCF ~{test_gvcf} differs in content from ~{truth_gvcf} by $DIFF_LINES lines"

      DIFF_LINES=$(diff <( gunzip -c -f ~{test_gvcf} | grep -v '^##' | cut -f 1-5,7- ) <( gunzip -c -f ~{truth_gvcf} | grep -v '^##' | cut -f 1-5,7- ) | grep -e "^<" | wc -l )

      echo "$DIFF_LINES" > diff_lines_quality.txt

      if [ $DIFF_LINES -eq 0 ]; then
        echo "However they ONLY differ in the quality column" >&2
      fi

    fi

    exit $exit_code
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 70 HDD"
    memory: "75 GiB"
    preemptible: 3
  }

  output {
    File gvcf_diff = "gvcf_diff.txt"
    File diff_lines_total = "diff_lines_total.txt"
    File diff_lines_quality = "diff_lines_quality.txt"
  }
}
