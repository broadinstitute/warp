version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifyRNAWithUMIs {

  input {
    Array[File] test_metrics
    Array[File] truth_metrics
    Array[File] test_text_metrics
    Array[File] truth_text_metrics
    File test_output_bam
    File truth_output_bam
    File test_transcriptome_bam
    File truth_transcriptome_bam
    File test_gene_tpm
    File truth_gene_tpm
    File test_gene_counts
    File truth_gene_counts
    File test_exon_counts
    File truth_exon_counts
    Boolean? done
  }

  call MetricsVerification.VerifyMetrics as CompareMetrics {
    input:
      test_metrics = test_metrics,
      truth_metrics = truth_metrics
  }

  call VerifyTasks.CompareTextFiles as CompareTextMetrics {
    input:
      test_text_files = test_text_metrics,
      truth_text_files = truth_text_metrics
  }

  call VerifyTasks.CompareBams as CompareOutputBam {
    input:
      test_bam = test_output_bam,
      truth_bam = truth_output_bam,
      lenient_header = true
  }

  call VerifyTasks.CompareTranscriptomeBams as CompareTranscriptomeBam {
    input:
      test_bam = test_transcriptome_bam,
      truth_bam = truth_transcriptome_bam,
      lenient_header = true
  }

  call VerifyTasks.CompareCompressedTextFiles as CompareGeneTpms {
    input:
      test_zip = test_gene_tpm,
      truth_zip = truth_gene_tpm
  }

  call VerifyTasks.CompareCompressedTextFiles as CompareGeneCounts {
    input:
      test_zip = test_gene_counts,
      truth_zip = truth_gene_counts
  }

  call VerifyTasks.CompareCompressedTextFiles as CompareExonCounts {
    input:
      test_zip = test_exon_counts,
      truth_zip = truth_exon_counts
  }

  output {
    Array[File] metric_comparison_report_files = CompareMetrics.metric_comparison_report_files
  }
  meta {
    allowNestedInputs: true
  }
}