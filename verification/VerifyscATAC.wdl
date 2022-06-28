version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyscATAC {
  input {
    File truth_bam
    File test_bam
    Array[File] truth_matrix_files
    Array[File] test_matrix_files
    
    Boolean? done
  }

  call VerifyTasks.CompareBams as CompareBams {
    input:
      test_bam       = test_bam,
      truth_bam      = truth_bam,
      lenient_header = true
  }

	call VerifyTasks.CompareTextFiles as CompareTextFiles{
		input:
			test_text_files  = test_matrix_files,
			truth_text_files = truth_matrix_files
	}

  output {}
}