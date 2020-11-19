version 1.0

import "../verification/VerifyTasks.wdl" as Tasks

## Copyright Broad Institute, 2018
##
## This WDL script is designed to verify (compare) the outputs of an MultiSampleArrays wdl.
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

workflow VerifyMultiSampleArrays {
  input {
    File truth_vcf
    File test_vcf

    Boolean fail_fast = true
  }

  call Tasks.CompareVcfs {
    input:
      file1 = truth_vcf,
      file2 = test_vcf,
      fail_fast = true
  }
  call SummarizeResults {
    input:
      compare_vcfs_exit_code = CompareVcfs.exit_code,
      compare_vcfs_results_file = CompareVcfs.report_file,
      fail_fast = true
  }
  output {
    Int exit_code = SummarizeResults.exit_code
    File report_file = SummarizeResults.report_file
  }
  meta {
    allowNestedInputs: true
  }
}

task SummarizeResults {
  input {
    Int compare_vcfs_exit_code
    File compare_vcfs_results_file
    Boolean fail_fast = true
  }

  command {
    fail_fast_value=~{true="true" false="" fail_fast}
    exit_code=~{compare_vcfs_exit_code}

    echo "Results of VerifyMultiSampleArrays Workflow:" >&2
    if [ "$exit_code" -eq "0" ]; then echo "Pass" >&2; else echo "Fail" >&2; fi

    echo >&2
    cat ~{compare_vcfs_results_file} >&2

    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stderr()
  }
}
