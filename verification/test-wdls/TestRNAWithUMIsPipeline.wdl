version 1.0

import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../verification/VerifyRNAWithUMIs.wdl" as VerifyRNAWithUMIs
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy
import "../../pipelines/broad/rna_seq/RNAWithUMIsPipeline.wdl" as RNAWithUMIsPipeline

workflow TestRNAWithUMIsPipeline {
  
  input {
      File? bam
      File? r1_fastq
      File? r2_fastq
      String read1Structure
      String read2Structure
      String output_basename
  
      String platform
      String library_name
      String platform_unit
      String read_group_name
      String sequencing_center = "BI"
  
      File starIndex
      File gtf
  
      File ref
      File refIndex
      File refDict
      File refFlat
      File ribosomalIntervals
      File exonBedFile
  
      File population_vcf
      File population_vcf_index
  
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

  call RNAWithUMIsPipeline.RNAWithUMIsPipeline {
    input:
        bam                  = bam,
        r1_fastq             = r1_fastq,
        r2_fastq             = r2_fastq,
        read1Structure       = read1Structure,
        read2Structure       = read2Structure,
        output_basename      = output_basename,
        platform             = platform,
        library_name         = library_name,
        platform_unit        = platform_unit,
        read_group_name      = read_group_name,
        sequencing_center    = sequencing_center,
        starIndex            = starIndex,
        gtf                  = gtf,
        ref                  = ref,
        refIndex             = refIndex,
        refDict              = refDict,
        refFlat              = refFlat,
        ribosomalIntervals   = ribosomalIntervals,
        exonBedFile          = exonBedFile,
        population_vcf       = population_vcf,
        population_vcf_index = population_vcf_index,
  }

  Array[String] pipeline_outputs = select_all([
                                    RNAWithUMIsPipeline.transcriptome_bam,
                                    RNAWithUMIsPipeline.output_bam,
                                    RNAWithUMIsPipeline.output_bam_index,
                                    RNAWithUMIsPipeline.rnaseqc2_gene_tpm,
                                    RNAWithUMIsPipeline.rnaseqc2_gene_counts,
                                    RNAWithUMIsPipeline.rnaseqc2_exon_counts,
                                    RNAWithUMIsPipeline.rnaseqc2_fragment_size_histogram,
                                    RNAWithUMIsPipeline.picard_insert_size_histogram,
                                    RNAWithUMIsPipeline.picard_base_distribution_by_cycle_pdf,
                                    RNAWithUMIsPipeline.picard_quality_by_cycle_pdf,
                                    RNAWithUMIsPipeline.picard_quality_distribution_pdf,
                                    RNAWithUMIsPipeline.fastqc_html_report
  ])

  Array[String] pipeline_metrics = select_all([
                                    RNAWithUMIsPipeline.transcriptome_duplicate_metrics,
                                    RNAWithUMIsPipeline.duplicate_metrics,
                                    RNAWithUMIsPipeline.picard_rna_metrics,
                                    RNAWithUMIsPipeline.picard_alignment_summary_metrics,
                                    RNAWithUMIsPipeline.picard_insert_size_metrics,
                                    RNAWithUMIsPipeline.picard_base_distribution_by_cycle_metrics,
                                    RNAWithUMIsPipeline.picard_quality_by_cycle_metrics,
                                    RNAWithUMIsPipeline.picard_quality_distribution_metrics
  ])

  Array[String] pipeline_text_metrics = select_all([RNAWithUMIsPipeline.rnaseqc2_metrics])

  #Copy results of pipeline to test results bucket
  call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics, pipeline_text_metrics]),
      vault_token_path          = vault_token_path,
      google_account_vault_path = google_account_vault_path,
      destination_cloud_path    = results_path
  }

  # If updating truth then copy pipeline results to truth bucket
  if (update_truth) {
    call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics, pipeline_text_metrics]),
      vault_token_path          = vault_token_path,
      google_account_vault_path = google_account_vault_path,
      destination_cloud_path    = truth_path
    }
  }

  if (!update_truth) {
    call Utilities.GetValidationInputs as GetMetricsInputs {
      input:
        input_files   = pipeline_metrics,
        results_path  = results_path,
        truth_path    = truth_path
    }

    call Utilities.GetValidationInputs as GetTextMetricsInputs {
      input:
        input_files   = pipeline_text_metrics,
        results_path  = results_path,
        truth_path    = truth_path
    }

    call Utilities.GetValidationInputs as GetBam {
      input:
        input_file    = RNAWithUMIsPipeline.output_bam,
        results_path  = results_path,
        truth_path    = truth_path
    }

    call Utilities.GetValidationInputs as GetTranscriptomeBam {
      input:
        input_file    = RNAWithUMIsPipeline.transcriptome_bam,
        results_path  = results_path,
        truth_path    = truth_path
    }

    call Utilities.GetValidationInputs as GetGeneTpm {
      input:
        input_file   = RNAWithUMIsPipeline.rnaseqc2_gene_tpm,
        results_path  = results_path,
        truth_path    = truth_path
    }

    call Utilities.GetValidationInputs as GetGeneCounts {
      input:
        input_file    = RNAWithUMIsPipeline.rnaseqc2_gene_counts,
        results_path  = results_path,
        truth_path    = truth_path
    }

    call Utilities.GetValidationInputs as GetExonCounts {
      input:
        input_file    = RNAWithUMIsPipeline.rnaseqc2_exon_counts,
        results_path  = results_path,
        truth_path    = truth_path
    }

    call VerifyRNAWithUMIs.VerifyRNAWithUMIs as Verify {
      input:
        test_metrics              = GetMetricsInputs.results_files,
        test_text_metrics         = GetTextMetricsInputs.results_files,
        test_output_bam           = GetBam.results_file,
        test_transcriptome_bam    = GetTranscriptomeBam.results_file,
        test_gene_tpm             = GetGeneTpm.results_file,
        test_gene_counts          = GetGeneCounts.results_file,
        test_exon_counts          = GetExonCounts.results_file,
        truth_metrics             = GetMetricsInputs.truth_files,
        truth_text_metrics        = GetTextMetricsInputs.truth_files,
        truth_output_bam          = GetBam.truth_file,
        truth_transcriptome_bam   = GetTranscriptomeBam.truth_file,
        truth_gene_tpm            = GetGeneTpm.truth_file,
        truth_gene_counts         = GetGeneCounts.truth_file,
        truth_exon_counts         = GetExonCounts.truth_file,
        done                      = CopyToTestResults.done
    }
  }

  output{}
}