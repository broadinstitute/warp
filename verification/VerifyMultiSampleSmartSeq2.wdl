version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyMultiSampleSmartSeq2 {
  input {
    Array[File] truth_bams
    Array[File] test_bams
    File truth_loom
    File test_loom

    Boolean? done
  }

  scatter (idx in range(length(truth_bams))){
    call VerifyTasks.CompareBams as CompareBams {
      input:
        test_bam       = test_bams[idx],
        truth_bam      = truth_bams[idx],
        lenient_header = true
    }
	}
  
  call VerifyTasks.CompareLooms as CompareLooms {
    input:
      test_loom  = test_loom,
      truth_loom = truth_loom
  }

  output{}
}