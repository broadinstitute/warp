version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyReblockGVCF {

  input {
    File test_gvcf
    File truth_gvcf

    Boolean fail_fast = true
  }
  
  call VerifyTasks.CompareVcfs {
    input:
      file1 = test_gvcf,
      file2 = truth_gvcf,
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

    echo "Results of VerifyReblockGVCF Workflow:" >&2
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