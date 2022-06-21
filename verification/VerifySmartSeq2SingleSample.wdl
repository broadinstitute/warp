version 1.0

import "../verification/VerifyMetrics.wdl" as VerifyMetrics
import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifySmartSeq2SingleSample {

  input {
    File truth_bam
    File test_bam
    File truth_transcriptome_bam
    File test_transcriptome_bam
    File truth_loom
    File test_loom
    Array[File] truth_metrics
    Array[File] test_metrics

    Boolean? done
  }
	call VerifyTasks.CompareBams as CompareBams {
    input:
      test_bam       = test_bam,
      truth_bam      = truth_bam,
      lenient_header = true
  }

	call VerifyTasks.CompareBams as CompareTranscriptomeBams {
    input:
      test_bam       = test_transcriptome_bam,
      truth_bam      = truth_transcriptome_bam,
      lenient_header = true
  }

	call VerifyMetrics.VerifyMetrics as CompareMetrics {
		input:
			test_metrics  = test_metrics,
			truth_metrics = truth_metrics
	}

	call VerifyTasks.CompareLooms as CompareLooms{
    input:
      test_loom  = test_loom,
      truth_loom = truth_loom
  }

    output {}
}