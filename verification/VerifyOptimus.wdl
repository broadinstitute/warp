version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyOptimus {
  
  input {
    File test_h5ad
    File truth_h5ad

    File test_bam
    File truth_bam

    File test_gene_metrics
    File truth_gene_metrics
    
    File test_cell_metrics
    File truth_cell_metrics

    File test_library_metrics
    File truth_library_metrics

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
  
  call VerifyTasks.CompareH5adFilesGEX as CompareH5adFilesOptimus {
        input:
            test_h5ad  = test_h5ad,
            truth_h5ad = truth_h5ad
    }

  call VerifyTasks.CompareLibraryFiles as CompareLibraryMetrics {
          input:
            test_text_file = test_library_metrics,
            truth_text_file = truth_library_metrics
  }
}
}