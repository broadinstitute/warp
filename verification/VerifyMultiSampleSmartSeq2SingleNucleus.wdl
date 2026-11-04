version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyMultiSampleSmartSeq2SingleNucleus {
	input {
    Array[File] truth_bams
    Array[File] test_bams
    File truth_h5ad
    File test_h5ad

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
  
  call VerifyTasks.CompareH5adFilesGEX as CompareH5adFiles {
    input:
      test_h5ad  = test_h5ad,
      truth_h5ad = truth_h5ad
  }

  output{}
}