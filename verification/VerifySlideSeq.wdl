version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifySlideSeq {

  input {
    File test_loom
    File truth_loom

    File test_bam
    File truth_bam

    File test_gene_metrics
    File truth_gene_metrics

    File test_cell_metrics
    File truth_cell_metrics

    Boolean? done
  }

  call VerifyTasks.CompareBams as CompareBams {
    input:
      test_bam       = test_bam,
      truth_bam      = truth_bam,
      lenient_header = true
  }

  call VerifyTasks.CompareCompressedTextFiles as CompareGeneMetrics {
    input:
      test_zip  = test_gene_metrics,
      truth_zip = truth_gene_metrics
  }

  call VerifyTasks.CompareCompressedTextFiles as CompareCellMetrics {
    input:
      test_zip  = test_cell_metrics,
      truth_zip = truth_cell_metrics
  }

  call VerifyTasks.CompareLooms as CompareLooms{
    input:
      test_loom  = test_loom,
      truth_loom = truth_loom
  }

}