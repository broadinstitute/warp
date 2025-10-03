version 1.0


import "../../pipelines/wdl/internal/rna_seq/BroadInternalRNAWithUMIs.wdl" as BroadInternalRNAWithUMIs
import "../../verification/VerifyRNAWithUMIs.wdl" as VerifyRNAWithUMIs
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestBroadInternalRNAWithUMIs {

    input {
      String reference_build
      String sample_lsid
      File r1_fastq
      File r2_fastq
      String read1Structure
      String read2Structure
      String output_basename
      String platform
      String library_name
      String platform_unit
      String read_group_name
      String sequencing_center = "BI"
      String? tdr_dataset_uuid
      String? tdr_sample_id

      # if there are very few duplicates, then relative change to duplication metrics can be high (0 vs 1), and some
      # metrics can be null (ESTIMATED_LIBRARY_SIZE if 0 duplicates, for example).  In these cases, just don't compare
      # transcriptome duplicate metrics
      Boolean compare_transcriptome_dup_metrics = true

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String environment
      String vault_token_path
      String vault_token_path_arrays
      String google_account_vault_path
    }

    meta {
      allowNestedInputs: true
    }
  
    call BroadInternalRNAWithUMIs.BroadInternalRNAWithUMIs {
      input:
        reference_build = reference_build,
        sample_lsid = sample_lsid,
        r1_fastq = r1_fastq,
        r2_fastq = r2_fastq,
        read1Structure = read1Structure,
        read2Structure = read2Structure,
        output_basename = output_basename,
        platform = platform,
        library_name = library_name,
        platform_unit = platform_unit,
        read_group_name = read_group_name,
        sequencing_center = sequencing_center,
        tdr_dataset_uuid = tdr_dataset_uuid,
        tdr_sample_id = tdr_sample_id,
        environment = environment,
        vault_token_path = vault_token_path_arrays
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    BroadInternalRNAWithUMIs.fastqc_html_report,
                                    BroadInternalRNAWithUMIs.picard_quality_distribution_pdf,
                                    BroadInternalRNAWithUMIs.picard_quality_by_cycle_pdf,
                                    BroadInternalRNAWithUMIs.picard_base_distribution_by_cycle_pdf,
                                    BroadInternalRNAWithUMIs.picard_insert_size_histogram,
                                    BroadInternalRNAWithUMIs.rnaseqc2_fragment_size_histogram,
                                    BroadInternalRNAWithUMIs.rnaseqc2_exon_counts,
                                    BroadInternalRNAWithUMIs.rnaseqc2_gene_counts,
                                    BroadInternalRNAWithUMIs.rnaseqc2_gene_tpm,
                                    BroadInternalRNAWithUMIs.output_bam_index,
                                    BroadInternalRNAWithUMIs.output_bam,
                                    BroadInternalRNAWithUMIs.transcriptome_bam,
                                    BroadInternalRNAWithUMIs.transcriptome_duplicate_metrics
                                    ],
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    BroadInternalRNAWithUMIs.picard_quality_distribution_metrics,
                                    BroadInternalRNAWithUMIs.picard_quality_by_cycle_metrics,
                                    BroadInternalRNAWithUMIs.picard_base_distribution_by_cycle_metrics,
                                    BroadInternalRNAWithUMIs.picard_insert_size_metrics,
                                    BroadInternalRNAWithUMIs.picard_alignment_summary_metrics,
                                    BroadInternalRNAWithUMIs.picard_rna_metrics,
                                    BroadInternalRNAWithUMIs.duplicate_metrics
                                    ],
                                    # File? outputs
                                    select_all([BroadInternalRNAWithUMIs.picard_fingerprint_detail_metrics]),
                                    select_all([BroadInternalRNAWithUMIs.picard_fingerprint_summary_metrics]),
                                    
    ])

    Array[String] pipeline_text_metrics = select_all([BroadInternalRNAWithUMIs.rnaseqc2_metrics, BroadInternalRNAWithUMIs.unified_metrics])

    # Copy results of pipeline to test results bucket
    call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics, pipeline_text_metrics]),
        vault_token_path          = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics, pipeline_text_metrics]),
          vault_token_path          = vault_token_path,
          google_account_vault_path = google_account_vault_path,
          destination_cloud_path    = truth_path
      }
    }

    # If not updating truth then we need to collect all input for the validation WDL
    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetTextMetrics {
          input:
            input_files = pipeline_text_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetOutputBam {
          input:
            input_file = BroadInternalRNAWithUMIs.output_bam,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetTranscriptomeBam {
          input:
            input_file = BroadInternalRNAWithUMIs.transcriptome_bam,
            results_path = results_path,
            truth_path = truth_path
        }
      call Utilities.GetValidationInputs as GetTranscriptomeDuplicationMetrics {
        input:
          input_file = BroadInternalRNAWithUMIs.transcriptome_duplicate_metrics,
          results_path  = results_path,
          truth_path    = truth_path
      }
        call Utilities.GetValidationInputs as GetGeneTpm {
          input:
            input_file = BroadInternalRNAWithUMIs.rnaseqc2_gene_tpm,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGeneCounts {
          input:
            input_file = BroadInternalRNAWithUMIs.rnaseqc2_gene_tpm,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetExonCounts {
          input:
            input_file = BroadInternalRNAWithUMIs.rnaseqc2_exon_counts,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyRNAWithUMIs.VerifyRNAWithUMIs as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_text_metrics = GetTextMetrics.truth_files,
          test_text_metrics = GetTextMetrics.results_files,
          truth_output_bam = GetOutputBam.truth_file, 
          test_output_bam = GetOutputBam.results_file,
          truth_transcriptome_bam = GetTranscriptomeBam.truth_file, 
          test_transcriptome_bam = GetTranscriptomeBam.results_file,
          test_transcriptome_duplicate_metrics = GetTranscriptomeDuplicationMetrics.results_file,
          truth_transcriptome_duplicate_metrics = GetTranscriptomeDuplicationMetrics.truth_file,
          truth_gene_tpm = GetGeneTpm.truth_file, 
          test_gene_tpm = GetGeneTpm.results_file,
          truth_gene_counts = GetGeneCounts.truth_file, 
          test_gene_counts = GetGeneCounts.results_file,
          truth_exon_counts = GetExonCounts.truth_file,
          test_exon_counts = GetExonCounts.results_file,
          transcriptome_deterministic = false,
          compare_transcriptome_dup_metrics = compare_transcriptome_dup_metrics,
          done = CopyToTestResults.done
      }
  
    }





}