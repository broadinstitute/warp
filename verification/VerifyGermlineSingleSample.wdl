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

    Boolean fail_fast = true
  }

  call MetricsVerification.VerifyMetrics as CompareMetrics {
    input:
      test_metrics = test_metrics,
      truth_metrics = truth_metrics,
      fail_fast = false
  }

  call Tasks.CompareGvcfs {
    input:
      test_gvcf = test_gvcf,
      truth_gvcf = truth_gvcf,
      fail_fast = false
  }

  call Tasks.CompareCrais {
    input:
      test_crai = test_crai,
      truth_crai = truth_crai,
      fail_fast = false
  }

  call Tasks.CompareCrams {
    input:
      test_cram = test_cram,
      test_crai = test_crai,
      truth_cram = truth_cram,
      truth_crai = truth_crai,
      fail_fast = false
  }

  call SummarizeResults {
    input:
      compare_gvcfs_exit_code = CompareGvcfs.exit_code,
      compare_gvcfs_results_file = CompareGvcfs.report_file,
      compare_crams_exit_code = CompareCrams.exit_code,
      compare_crams_results_file = CompareCrams.report_file,
      compare_crais_exit_code = CompareCrais.exit_code,
      compare_crais_results_file = CompareCrais.report_file,
      compare_metrics_exit_code = CompareMetrics.exit_code,
      compare_metrics_results_file = CompareMetrics.report_file,
      fail_fast = true
  }
  output {
    Int exit_code = SummarizeResults.exit_code
    File report_file = SummarizeResults.report_file
    Array[File] metric_comparison_report_files = CompareMetrics.metric_comparison_report_files
  }
  meta {
    allowNestedInputs: true
  }
}

task SummarizeResults {
  input {
    Int compare_gvcfs_exit_code
    File compare_gvcfs_results_file
    Int compare_crams_exit_code
    File compare_crams_results_file
    Int compare_crais_exit_code
    File compare_crais_results_file
    Int compare_metrics_exit_code
    File compare_metrics_results_file
    Boolean fail_fast = true
  }

  command {
    fail_fast_value=~{true="true" false="" fail_fast}
    exit_code=$((~{compare_gvcfs_exit_code} + ~{compare_crams_exit_code} + ~{compare_crais_exit_code} + ~{compare_metrics_exit_code}))

    echo "Results of VerifyGermlineSingleSample Workflow:"
    if [ "$exit_code" -eq "0" ]; then echo "Pass"; else echo "Fail"; fi

    echo
    cat ~{compare_gvcfs_results_file}
    echo
    cat ~{compare_crams_results_file}
    echo
    cat ~{compare_crais_results_file}
    echo
    cat ~{compare_metrics_results_file}

    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 70 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stdout()
  }
}
