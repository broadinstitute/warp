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
        file1 = truth_metrics[idx],
        file2 = test_metrics[idx],
        output_file = "metric_~{idx}.txt",
        metrics_to_ignore = []
    }
  }

  # call ConsolidateErrors as ConsolidateErrors {
  #   input:
  #     error_files = CompareMetricFiles.report_file   
  # }

  output {
    Array[File] metric_comparison_report_files = CompareMetricFiles.report_file
    #Consolidate failed_metrics to just one file:
    File failed_metrics_file = write_tsv(CompareMetricFiles.failed_metrics)
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
      ~{true="--METRICS_TO_IGNORE" false="" length(metrics_to_ignore) > 0} ~{default="" sep=" --METRICS_TO_IGNORE " metrics_to_ignore};

      #TODO: Check if file has the text "Metrics are NOT equal", if it does: create a file that contains just the first 3 lines fo the output file
      #If the text is not found, ie. metrics are equal, do nothing or write an empty file.

      # Check for the string "Metrics are NOT equal"
      if grep -q "Metrics are NOT equal" ~{output_file}
      then
          # If string exists, copy output_file to failed_metrics_file.txt
          cat ~{output_file} > failed_metrics_file.txt
      else
          # If string does not exist, create an empty file named failed_metrics_file.txt
          touch failed_metrics_file.txt
      fi
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
  }
  output {
    File report_file = output_file
    Array[String] failed_metrics = read_lines("failed_metrics_file.txt")
  }
}
