version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as Tasks
import "../verification/VerifyNA12878.wdl" as VerifyNA12878

workflow VerifyUltimaGenomicsWholeGenomeGermline {

  input {
    Array[File] test_metrics
    Array[File] truth_metrics

    File test_cram
    File truth_cram
    File test_crai
    File truth_crai

    File test_vcf
    File truth_vcf
    File test_filtered_vcf
    File test_filtered_vcf_index
    File truth_filtered_vcf
    File truth_filtered_vcf_index

    File test_gvcf
    File truth_gvcf
    File truth_gvcf_index
    File test_gvcf_index

    String sample_name

    Boolean? done
  }

  ## TEMP - Can't use this simple call, because the picard version used by UltimaGenomicsWholeGenomeGermline is based on an old version of Picard,
#  call MetricsVerification.VerifyMetrics as CompareMetrics {
#    input:
#      test_metrics = test_metrics,
#      truth_metrics = truth_metrics
#  }
  # So, I have to do all this - when the UltimaGenomicsWholeGenomeGermline picard is properly rebased, this should all be able to be rolled back.
  call MetricsVerification.CompareTwoNumbers {
    input:
      num1 = length(test_metrics),
      num2 = length(truth_metrics),
      error_msg = "Different number of metric files"
  }

  String alignment_metrics_ext = ".alignment_summary_metrics"
  String duplicate_metrics_ext = ".duplicate_metrics"

  scatter (idx in range(length(truth_metrics))) {
    String metrics_basename = basename(truth_metrics[idx])
    Boolean is_alignment_metrics_file = basename(metrics_basename, alignment_metrics_ext) != metrics_basename
    Boolean is_duplicate_metrics_file = basename(metrics_basename, duplicate_metrics_ext) != metrics_basename
    if ((!is_alignment_metrics_file) && (!is_duplicate_metrics_file)) {
      call MetricsVerification.CompareMetricFiles {
        input:
          dependency_input = CompareTwoNumbers.output_file,
          file1 = test_metrics[idx],
          file2 = truth_metrics[idx],
          output_file = "metric_~{idx}.txt",
          metrics_to_ignore = []
      }
    }
    if ((is_alignment_metrics_file) || (is_duplicate_metrics_file)) {
      call CompareOldMetricFiles {
        input:
          dependency_input = CompareTwoNumbers.output_file,
          file1 = test_metrics[idx],
          file2 = truth_metrics[idx],
          output_file = "metric_~{idx}.txt",
          metrics_to_ignore = []
      }
    }
  }
  ## END TEMP

  call Tasks.CompareCrams {
    input:
      test_cram = test_cram,
      test_crai = test_crai,
      truth_cram = truth_cram,
      truth_crai = truth_crai
  }

  call Tasks.CompareCrais {
    input:
      test_crai = test_crai,
      truth_crai = truth_crai
  }

  call Tasks.CompareVcfs as CompareVcfs {
    input:
      file1 = test_vcf,
      file2 = truth_vcf,
      patternForLinesToExcludeFromComparison = "^##GATKCommandLine"
  }

  call Tasks.CompareVcfs as CompareFilteredVcfs {
    input:
      file1 = test_filtered_vcf,
      file2 = truth_filtered_vcf,
      patternForLinesToExcludeFromComparison = "^##"
  }

  call Tasks.CompareVcfs as CompareGvcfs {
    input:
      file1 = test_gvcf,
      file2 = truth_gvcf,
      patternForLinesToExcludeFromComparison = "^##"
  }

  call Tasks.CompareVCFsVerbosely as CompareGVcfsVerbosely {
    input:
      actual = test_gvcf,
      expected = truth_gvcf,
      actual_index = test_gvcf_index,
      expected_index = truth_gvcf_index,
      extra_args = " --ignore-attribute AVERAGE_TREE_SCORE --ignore-attribute CALIBRATION_SENSITIVITY "
                   + "--ignore-attribute TREE_SCORE --ignore-attribute SCORE --ignore-filters --ignore-attribute AS_FilterStatus "
                     + "--ignore-attribute ExcessHet --ignore-star-attributes --allow-nan-mismatch --ignore-attribute END"
  }

  call VerifyNA12878.VerifyNA12878 {
    input:
      vcf_files = [test_filtered_vcf, truth_filtered_vcf],
      vcf_file_indexes = [test_filtered_vcf_index, truth_filtered_vcf_index],
      vcf_names = ["test","truth"],
      sample_name = sample_name
  }

  meta {
    allowNestedInputs: true
  }

  output{}
}

task CompareGvcfs {

  input {
    File test_gvcf
    File truth_gvcf

  }

  command {
    exit_code=0

    DIFF_LINES=$(diff <(gunzip -c -f ~{test_gvcf} | grep -v '^##') <(gunzip -c -f ~{truth_gvcf} | grep -v '^##') | grep -e "^<" | wc -l)
    if [ $DIFF_LINES -ge 10 ]; then
      exit_code=1
      echo "Error: GVCF ~{test_gvcf} differs in content from ~{truth_gvcf} by $DIFF_LINES lines" >&2
      DIFF_LINES=$(diff <(gunzip -c -f ~{test_gvcf} | grep -v '^##' | cut -f 1-5,7-) <(gunzip -c -f ~{truth_gvcf} | grep -v '^##' | cut -f 1-5,7-) | grep -e "^<" | wc -l)
      if [ $DIFF_LINES -eq 0 ]; then
        echo "However they ONLY differ in the quality column" >&2
      fi
    fi
    exit $exit_code
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 70 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  
}

task CompareOldMetricFiles {
  input {
    File? dependency_input
    File file1
    File file2
    String output_file
    Array[String] metrics_to_ignore
  }

  command <<<
    java -Xmx3g -Dpicard.useLegacyParser=false  -jar /usr/gitc/picard.jar \
    CompareMetrics \
    --INPUT ~{file1} \
    --INPUT ~{file2} \
    --OUTPUT ~{output_file} \
    ~{true="--METRICS_TO_IGNORE" false="" length(metrics_to_ignore) > 0} ~{default="" sep=" --METRICS_TO_IGNORE " metrics_to_ignore}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:samtools_picard_bwa_snapshot_UG"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
  }
  output {
    File report_file = output_file
  }
}
