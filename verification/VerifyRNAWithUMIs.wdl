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
    File test_transcriptome_duplicate_metrics
    File truth_transcriptome_duplicate_metrics
    Boolean transcriptome_deterministic
    Boolean compare_transcriptome_dup_metrics

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

  if (transcriptome_deterministic) {
    call VerifyTasks.CompareBams as CompareTranscriptomeBam {
      input:
        test_bam = test_transcriptome_bam,
        truth_bam = truth_transcriptome_bam,
        lenient_header = true
    }
  }

  if (!transcriptome_deterministic) {
    call RemoveMultiMappers as RemoveMultiMappersTestBam{
      input:
        bam = test_transcriptome_bam
    }

    call RemoveMultiMappers as RemoveMultiMappersTruthBam{
      input:
        bam = truth_transcriptome_bam
    }

    call CompareTranscriptomeBamNonDeterministic {
      input:
        test_bam = RemoveMultiMappersTestBam.output_bam,
        truth_bam = RemoveMultiMappersTruthBam.output_bam
    }

    call CheckTranscriptomeBamComparisonWithTolerance {
      input:
        comparison = CompareTranscriptomeBamNonDeterministic.comparison,
        tolerance = 0.006
    }
  }

  if (compare_transcriptome_dup_metrics) {
    call MetricsVerification.CompareMetricFiles as CompareTranscriptomeDuplicationMetrics {
      input:
        file1 = truth_transcriptome_duplicate_metrics,
        file2 = test_transcriptome_duplicate_metrics,
        output_file  = "transcriptome_duplication_metrics_comparison.txt",
        extra_args = if transcriptome_deterministic then [] else ["--METRIC_ALLOWABLE_RELATIVE_CHANGE READ_PAIR_DUPLICATES:0.005",
                                                                  "--METRIC_ALLOWABLE_RELATIVE_CHANGE READ_PAIR_OPTICAL_DUPLICATES:0.01",
                                                                  "--METRIC_ALLOWABLE_RELATIVE_CHANGE PERCENT_DUPLICATION:0.005",
                                                                  "--METRIC_ALLOWABLE_RELATIVE_CHANGE ESTIMATED_LIBRARY_SIZE:0.005",
                                                                  "--IGNORE_HISTOGRAM_DIFFERENCES"]
    }
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

task RemoveMultiMappers {
  input {
    File bam
  }

  String basename = basename(bam, ".bam")
  Int disk_size = 2*ceil(size(bam,"GB")) + 20

  command <<<
    set -e
    set -o pipefail

    samtools view -q 4 -bo ~{basename}.multimapper_removed.bam ~{bam}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_bam = "~{basename}.multimapper_removed.bam"
  }
}

task CompareTranscriptomeBamNonDeterministic {
  input {
    File test_bam
    File truth_bam
  }

  Float bam_size = size(test_bam, "GiB") + size(truth_bam, "GiB")
  Int disk_size = ceil(bam_size * 4) + 20

  command <<<
    set -e
    set -o pipefail

    java -Xms3500m -Xmx7000m -jar /usr/gitc/picard.jar \
    CompareSAMs \
    ~{test_bam} \
    ~{truth_bam} \
    O=comparison.tsv \
    LENIENT_HEADER=true \
    LENIENT_LOW_MQ_ALIGNMENT=true
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.27.5-1667410556"
    disks: "local-disk " + disk_size + " HDD"
    cpu: 2
    memory: "7500 MiB"
    preemptible: 3
    continueOnReturnCode: [0, 1]
  }

  output {
    File comparison = "comparison.tsv"
  }
}

task CheckTranscriptomeBamComparisonWithTolerance {
  input {
    File comparison
    Float tolerance
  }

  command <<<
    set -e
    set -o pipefail

    python3 << EOF
    import csv
    # Picard CompareMetrics output: 5 preamble lines, then header (row 5), then data (row 6).
    with open("~{comparison}") as fh:
        row = dict(zip(*list(csv.reader(fh, delimiter="\t"))[5:7]))

    frac_missing = (float(row['MISSING_LEFT']) + float(row['MISSING_RIGHT'])) / float(row['MAPPINGS_MATCH'])
    assert frac_missing < ~{tolerance}, f'frac missing is {frac_missing} which is greater than tolerance of ~{tolerance}'
    assert float(row['MAPPINGS_DIFFER'])==0, f"{row['MAPPINGS_DIFFER']} mappings differ"
    assert float(row['UNMAPPED_LEFT'])==0, f"{row['UNMAPPED_LEFT']} unmapped in left file"
    assert float(row['UNMAPPED_RIGHT'])==0, f"{row['UNMAPPED_RIGHT']} unmapped in right file"
    assert float(row['DUPLICATE_MARKINGS_DIFFER'])==0, f"{row['DUPLICATE_MARKINGS_DIFFER']} duplicate markings differ (all duplicates should have been removed)"

    EOF
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 3
  }
}