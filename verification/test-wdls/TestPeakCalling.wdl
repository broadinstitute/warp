version 1.0


import "../../pipelines/wdl/peak_calling/PeakCalling.wdl" as PeakCalling
import "../../verification/VerifyPeakCalling.wdl" as VerifyPeakCalling
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestPeakCalling {

    input {
      File annotations_gtf
      File metrics_h5ad
      File chrom_sizes
      String output_base_name
      Int min_counts = 5000
      Int min_tsse = 10
      Int max_counts = 100000
      Float probability_threshold = 0.5
      Int disk_size = 500
      Int mem_size = 64
      Int nthreads = 4
      String cloud_provider

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call PeakCalling.PeakCalling {
      input:
        annotations_gtf = annotations_gtf,
        metrics_h5ad = metrics_h5ad,
        chrom_sizes = chrom_sizes,
        output_base_name = output_base_name,
        min_counts = min_counts,
        min_tsse = min_tsse,
        max_counts = max_counts,
        probability_threshold = probability_threshold,
        disk_size = disk_size,
        mem_size = mem_size,
        nthreads = nthreads,
        cloud_provider = cloud_provider
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    PeakCalling.cellbypeak_h5ad,
                                    PeakCalling.cellbybin_h5ad,
                                    ],
                                    
    ])

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs]),
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs]),
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetCellbybinH5ad {
          input:
            input_file = PeakCalling.cellbybin_h5ad,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCellbypeakH5ad {
          input:
            input_file = PeakCalling.cellbypeak_h5ad,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyPeakCalling.VerifyPeakCalling as Verify {
        input:
          truth_cellbybin_h5ad = GetCellbybinH5ad.truth_file, 
          test_cellbybin_h5ad = GetCellbybinH5ad.results_file,
          truth_cellbypeak_h5ad = GetCellbypeakH5ad.truth_file, 
          test_cellbypeak_h5ad = GetCellbypeakH5ad.results_file,
          done = CopyToTestResults.done
      }
    }
}