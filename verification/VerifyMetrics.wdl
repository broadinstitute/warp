version 1.0

workflow VerifyMetrics {

  input {
    Array[File] test_metrics
    Array[File] truth_metrics
  }

  call CompareTwoNumbers as CompareNumberOfMetricFiles {
    input:
      num1 = length(test_metrics),
      num2 = length(truth_metrics),
      error_msg = "Different number of metric files"
  }

  scatter (idx in range(length(truth_metrics))) {
    call CompareMetricFiles {
      input:
        dependency_input = CompareNumberOfMetricFiles.output_file,
        file1 = test_metrics[idx],
        file2 = truth_metrics[idx],
        output_file = "metric_~{idx}.txt",
        metrics_to_ignore = []
    }
  }

  output {
    Array[File] metric_comparison_report_files = CompareMetricFiles.report_file
  }
  meta {
    allowNestedInputs: true
  }
}

task CompareTwoNumbers {
  input {
    Int num1
    Int num2
    String? error_msg = "Numbers differ"
  }

  command {
    echo "Numbers agree" > output_file.txt
    if [ ~{num1} -ne ~{num2} ]; then
        echo "Error: ~{error_msg} (~{num1} vs. ~{num2})." > output_file.txt
        exit 1
    fi
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
  }
  output {
    File output_file = "output_file.txt"
  }
}

task CompareMetricFiles {
  input {
    File? dependency_input
    File file1
    File file2
    String output_file
    Array[String] metrics_to_ignore
  }

  command <<<
    java -Xmx3g -Dpicard.useLegacyParser=false  -jar /usr/picard/picard.jar \
      CompareMetrics \
      --INPUT ~{file1} \
      --INPUT ~{file2} \
      --OUTPUT ~{output_file} \
      ~{true="--METRICS_TO_IGNORE" false="" length(metrics_to_ignore) > 0} ~{default="" sep=" --METRICS_TO_IGNORE " metrics_to_ignore}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
  }
  output {
    File report_file = output_file
  }
}