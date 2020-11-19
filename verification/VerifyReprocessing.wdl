version 1.0

import "../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample

workflow VerifyReprocessing {

  input {
    Array[BamPair] bam_pairs

    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_cram
    File truth_crai
    File test_cram
    File test_crai

    File truth_gvcf
    File test_gvcf
  }

  scatter(pair in bam_pairs) {
    call CompareBams {
      input:
        test_bam = pair.test_bam,
        truth_bam = pair.truth_bam,
        fail_fast = false
    }
  }

  call VerifyGermlineSingleSample.VerifyGermlineSingleSample {
    input:
      truth_metrics = truth_metrics,
      test_metrics = test_metrics,
      truth_cram = truth_cram,
      truth_crai = truth_crai,
      test_cram = test_cram,
      test_crai = test_crai,
      truth_gvcf = truth_gvcf,
      test_gvcf = test_gvcf,
      fail_fast = false
  }

  call SummarizeResults {
    input:
      compare_bams_exit_code = CompareBams.exit_code,
      compare_bams_results_file = CompareBams.report_file,
      verify_germline_single_sample_exit_code = VerifyGermlineSingleSample.exit_code,
      verify_germline_single_sample_results_file = VerifyGermlineSingleSample.report_file,
      fail_fast = true
  }
  output {
    Int exit_code = SummarizeResults.exit_code
    File report_file = SummarizeResults.report_file
    Array[File] metric_comparison_report_files = VerifyGermlineSingleSample.metric_comparison_report_files
  }
  meta {
    allowNestedInputs: true
  }
}

struct BamPair {
  File test_bam
  File truth_bam
}

task SummarizeResults {
  input {
    Array[Int] compare_bams_exit_code
    Array[File] compare_bams_results_file
    Int verify_germline_single_sample_exit_code
    File verify_germline_single_sample_results_file
    Boolean fail_fast = true
  }

  command <<<
    fail_fast_value=~{true="true" false="" fail_fast}
    mapfile -t files < ~{write_lines(compare_bams_results_file)}
    mapfile -t exit_codes < ~{write_lines(compare_bams_exit_code)}

    echo "Results of VerifyReprocessing Workflow:"

    exit_code=~{verify_germline_single_sample_exit_code}
    for ((i=0;i<${#files[@]};++i)); do
      echo "------------"
      cat ${files[$i]}
      if [ "${exit_codes[$i]}" -ne "0" ]; then
        exit_code=${exit_codes[$i]}
      fi
    done

    echo
    cat ~{verify_germline_single_sample_results_file}
    echo

    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  >>>
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stdout()
  }
}

task CompareBams {

  input {
    File test_bam
    File truth_bam
    Boolean fail_fast = true
  }

  Float bam_size = size(test_bam, "GiB") + size(truth_bam, "GiB")
  Int disk_size = ceil(bam_size * 4) + 20

  command {
    fail_fast_value=~{true="true" false="" fail_fast}

    echo "Results of CompareBams:"
    echo -e "Test:\t~{test_bam}"
    echo -e "Truth:\t~{truth_bam}"

    # ApplyBQSR changes quality scores, so we need to strip those out of both BAMs
    cmp \
      <(samtools sort -n ~{test_bam} | samtools view | cut -d$'\t' -f 1-10,12-) \
      <(samtools sort -n ~{truth_bam} | samtools view | cut -d$'\t' -f 1-10,12-) > cmp_out.txt

    exit_code=$?
    if [ "$exit_code" -eq "0" ]; then
      echo -e "Pass\tBAMs do not differ"
    else
      echo -e "Fail\tBAMs differ\t`cat cmp_out.txt`"
    fi
    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    disks: "local-disk " + disk_size + " HDD"
    cpu: 2
    memory: "7.5 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stdout()
  }
}
