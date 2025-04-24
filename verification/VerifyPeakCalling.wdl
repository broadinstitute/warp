version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyPeakCalling {

    input {
        File test_cellbybin_h5ad
        File truth_cellbybin_h5ad

        File test_cellbypeak_h5ad
        File truth_cellbypeak_h5ad

        Boolean? done
    }

    call VerifyTasks.CompareH5adFilesATAC as CompareH5adFilesCellByBin {
        input:
            test_h5ad  = test_cellbybin_h5ad,
            truth_h5ad = truth_cellbybin_h5ad
    }
    call VerifyTasks.CompareH5adFilesATAC as CompareH5adFilesCellByPeak {
        input:
            test_h5ad  = test_cellbypeak_h5ad,
            truth_h5ad = truth_cellbypeak_h5ad
    }
}