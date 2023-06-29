version 1.0


import "../../pipelines/skylab/smartseq2_single_sample/SmartSeq2SingleSample.wdl" as SmartSeq2SingleSample
import "../../verification/VerifySmartSeq2SingleSample.wdl" as VerifySmartSeq2SingleSample
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestSmartSeq2SingleSample {

    input {
      File genome_ref_fasta
      File rrna_intervals
      File gene_ref_flat
      File hisat2_ref_index
      File hisat2_ref_trans_index
      File rsem_ref_index
      String hisat2_ref_name
      String hisat2_ref_trans_name
      String stranded
      String input_id
      String? input_name
      String? input_id_metadata_field
      String? input_name_metadata_field
      String output_name
      File fastq1
      File? fastq2
      Boolean paired_end

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
  
    call SmartSeq2SingleSample.SmartSeq2SingleSample {
      input:
        genome_ref_fasta = genome_ref_fasta,
        rrna_intervals = rrna_intervals,
        gene_ref_flat = gene_ref_flat,
        hisat2_ref_index = hisat2_ref_index,
        hisat2_ref_trans_index = hisat2_ref_trans_index,
        rsem_ref_index = rsem_ref_index,
        hisat2_ref_name = hisat2_ref_name,
        hisat2_ref_trans_name = hisat2_ref_trans_name,
        stranded = stranded,
        input_id = input_id,
        input_name = input_name,
        input_id_metadata_field = input_id_metadata_field,
        input_name_metadata_field = input_name_metadata_field,
        output_name = output_name,
        fastq1 = fastq1,
        fastq2 = fastq2,
        paired_end = paired_end
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    SmartSeq2SingleSample.loom_output_files,
                                    SmartSeq2SingleSample.rsem_isoform_results,
                                    SmartSeq2SingleSample.rsem_gene_results,
                                    SmartSeq2SingleSample.aligned_transcriptome_bam,
                                    SmartSeq2SingleSample.bam_index,
                                    SmartSeq2SingleSample.aligned_bam,
                                    ],
                                    # Array[File] outputs
                                    SmartSeq2SingleSample.group_results,
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    SmartSeq2SingleSample.rna_metrics,
                                    SmartSeq2SingleSample.bait_bias_summary_metrics,
                                    SmartSeq2SingleSample.quality_by_cycle_metrics,
                                    SmartSeq2SingleSample.quality_distribution_metrics,
                                    ],
                                    # File? outputs
                                    select_all([SmartSeq2SingleSample.insert_size_metrics]),
                                    
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
        call Utilities.GetValidationInputs as GetBam {
          input:
            input_file = SmartSeq2SingleSample.aligned_bam,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetTranscriptomeBam {
          input:
            input_file = SmartSeq2SingleSample.aligned_transcriptome_bam,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetLoom {
          input:
            input_file = SmartSeq2SingleSample.loom_output_files,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifySmartSeq2SingleSample.VerifySmartSeq2SingleSample as Verify {
        input:
          truth_bam = GetBam.truth_file, 
          test_bam = GetBam.results_file,
          truth_transcriptome_bam = GetTranscriptomeBam.truth_file, 
          test_transcriptome_bam = GetTranscriptomeBam.results_file,
          truth_loom = GetLoom.truth_file, 
          test_loom = GetLoom.results_file,
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          done = CopyToTestResults.done
      }
    }
}