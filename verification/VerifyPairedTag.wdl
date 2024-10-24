version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyPairedTag {

    input {
        File test_optimus_h5ad
        File truth_optimus_h5ad

        File test_optimus_bam
        File truth_optimus_bam

        File test_gene_metrics
        File truth_gene_metrics

        File test_cell_metrics
        File truth_cell_metrics

        File test_atac_bam
        File truth_atac_bam

        File test_fragment_file
        File truth_fragment_file

        File test_atac_h5ad
        File truth_atac_h5ad

        File test_library_metrics
        File truth_library_metrics

        Boolean? done
    }

    call VerifyTasks.CompareBams as CompareOptimusBams {
        input:
            test_bam       = test_optimus_bam,
            truth_bam      = truth_optimus_bam,
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

    call VerifyTasks.CompareBams as CompareAtacBams {
        input:
            test_bam       = test_atac_bam,
            truth_bam      = truth_atac_bam,
            lenient_header = true
    }
    call VerifyTasks.CompareTabix as CompareFragment {
        input:
            test_fragment_file  = test_fragment_file,
            truth_fragment_file = truth_fragment_file          
    }
    call VerifyTasks.CompareH5adFilesATAC as CompareH5adFilesATAC {
        input:
            test_h5ad  = test_atac_h5ad,
            truth_h5ad = truth_atac_h5ad
    }
    call VerifyTasks.CompareH5adFilesGEX as CompareH5adFilesOptimus {
        input:
            test_h5ad  = test_optimus_h5ad,
            truth_h5ad = truth_optimus_h5ad
    }
    call VerifyTasks.CompareLibraryFiles as CompareLibraryMetrics {
        input:
            test_text_file = test_library_metrics,
            truth_text_file = truth_library_metrics
    }
}