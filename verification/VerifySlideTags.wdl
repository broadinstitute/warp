version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifySlideTags {

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

        # spatial and positioning outputs 
        File test_spatial_output_h5
        File truth_spatial_output_h5

        File test_seurat_qs 
        File truth_seurat_qs 

        File test_coords_csv
        File truth_coords_csv

        File test_coords2_csv
        File truth_coords2_csv

        File test_intermediates_file
        File truth_intermediates_file

        Boolean? done
  }
    
    call VerifyTasks.CompareBams as CompareOptimusBams {
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
    
    call VerifyTasks.CompareH5Files as CompareSpatialOutputH5 {
        input:
            test_h5  = test_spatial_output_h5,
            truth_h5 = truth_spatial_output_h5
    }

    call VerifyTasks.CompareTextFiles as CompareCSV {
	      input:
		        test_text_files  = [test_coords_csv],
		        truth_text_files = [truth_coords_csv]
    }
    
    call VerifyTasks.CompareTextFiles as CompareCSV2 {
	      input:
		        test_text_files  = [test_coords2_csv],
		        truth_text_files = [truth_coords2_csv]
    }

    call VerifyTasks.CompareCompressedTextFiles as CompareTAR {
	      input:
		        test_zip  = test_intermediates_file,
		        truth_zip = truth_intermediates_file
    }

}
