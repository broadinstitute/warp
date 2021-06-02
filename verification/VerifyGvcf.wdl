version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyGvcf {

  input {
    File test_gvcf
    File truth_gvcf
  }
  
  call VerifyTasks.CompareVcfs {
    input:
      file1 = test_gvcf,
      file2 = truth_gvcf
  }
  meta {
    allowNestedInputs: true
  }
}