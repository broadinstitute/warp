version 1.0


import "../../pipelines/skylab/slideseq/SlideSeq.wdl" as SlideSeq
import "../../verification/VerifySlideSeq.wdl" as VerifySlideSeq
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestSlideSeq {

    input {
      Array[File] r1_fastq
      Array[File] r2_fastq
      Array[File]? i1_fastq
      String input_id
      String read_structure
      File tar_star_reference
      File annotations_gtf
      String output_bam_basename
      Boolean count_exons
      File bead_locations

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String vault_token_path
      String google_account_vault_path
    }

    meta {
      allowNestedInputs: true
    }
  
    call SlideSeq.SlideSeq {
      input:
        r1_fastq = r1_fastq,
        r2_fastq = r2_fastq,
        i1_fastq = i1_fastq,
        input_id = input_id,
        read_structure = read_structure,
        tar_star_reference = tar_star_reference,
        annotations_gtf = annotations_gtf,
        output_bam_basename = output_bam_basename,
        count_exons = count_exons,
        bead_locations = bead_locations
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    SlideSeq.matrix_col_index,
                                    SlideSeq.matrix_row_index,
                                    SlideSeq.matrix,
                                    SlideSeq.bam,
                                    ],
                                    # File? outputs
                                    select_all([SlideSeq.loom_output_file]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    SlideSeq.gene_metrics,
                                    SlideSeq.cell_metrics,
                                    SlideSeq.umi_metrics,
                                    ],
                                    
    ])

    # Copy results of pipeline to test results bucket
    call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        vault_token_path          = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
          vault_token_path          = vault_token_path,
          google_account_vault_path = google_account_vault_path,
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
          call Utilities.GetValidationInputs as GetLoom {
            input:
              input_file = SlideSeq.loom_output_file,
              results_path = results_path,
              truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetBam {
          input:
            input_file = SlideSeq.bam,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGeneMetrics {
          input:
            input_file = SlideSeq.gene_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCellMetrics {
          input:
            input_file = SlideSeq.cell_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetUmiMetrics {
            input:
                input_file = SlideSeq.umi_metrics,
                results_path = results_path,
                truth_path = truth_path
        }

      call VerifySlideSeq.VerifySlideSeq as Verify {
        input:
          truth_loom = GetLoom.truth_file,
          test_loom = GetLoom.results_file,
          truth_bam = GetBam.truth_file, 
          test_bam = GetBam.results_file,
          truth_gene_metrics = GetGeneMetrics.truth_file, 
          test_gene_metrics = GetGeneMetrics.results_file,
          truth_cell_metrics = GetCellMetrics.truth_file, 
          test_cell_metrics = GetCellMetrics.results_file,
          truth_umi_metrics = GetUmiMetrics.truth_file,
          test_umi_metrics = GetUmiMetrics.results_file,
          done = CopyToTestResults.done
      }
    }
}