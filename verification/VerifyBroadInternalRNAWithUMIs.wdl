version 1.0

import "../verification/VerifyRNAWithUMIs.wdl" as VerifyRNAWithUMIs

workflow VerifyBroadInternalRNAWithUMIs {

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
  }

  call VerifyRNAWithUMIs.VerifyRNAWithUMIs {
    input:
      test_metrics = test_metrics,
      truth_metrics = truth_metrics,
      test_text_metrics = test_text_metrics,
      truth_text_metrics = truth_text_metrics,
      test_output_bam = test_output_bam,
      truth_output_bam = truth_output_bam,
      test_transcriptome_bam = test_transcriptome_bam,
      truth_transcriptome_bam = truth_transcriptome_bam,
      test_gene_tpm = test_gene_tpm,
      truth_gene_tpm = truth_gene_tpm,
      test_gene_counts = test_gene_counts,
      truth_gene_counts = truth_gene_counts,
      test_exon_counts = test_exon_counts,
      truth_exon_counts = truth_exon_counts
  }

  output {
    Array[File] metric_comparison_report_files = VerifyRNAWithUMIs.metric_comparison_report_files
  }
  meta {
    allowNestedInputs: true
  }
}