version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyGvcf {

  input {
    File test_gvcf
    File test_gvcf_index
    File truth_gvcf
    File truth_gvcf_index

    Boolean? done
  }

  call VerifyTasks.CompareVCFsVerbosely {
    input:
      actual = test_gvcf,
      actual_index = test_gvcf_index,
      expected = truth_gvcf,
      expected_index = truth_gvcf_index
  }

  call VerifyTasks.CompareVcfs {
    input:
      file1 = test_gvcf,
      file2 = truth_gvcf,
      patternForLinesToExcludeFromComparison = "^##GATKCommandLine"
  }
  meta {
    allowNestedInputs: true
  }
}