version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks
import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyGermlineSingleSample.wdl" as GermlineVerification
import "../verification/VerifyNA12878.wdl" as VerifyNA12878

workflow VerifyJointGenotyping {

  input {
    Array[File] test_vcfs
    Array[File] test_vcf_indexes
    Array[File] truth_vcfs
    Array[File] truth_vcf_indexes

    Array[File] test_intervals
    Array[File] truth_intervals

    Array[File] test_metrics
    Array[File] truth_metrics

    File test_fingerprint
    File truth_fingerprint

    Boolean? done
  }

  scatter (idx in range(length(truth_vcfs))) {
    call GermlineVerification.CompareGvcfs {
      input:
        test_gvcf = test_vcfs[idx],
        truth_gvcf = truth_vcfs[idx]
    }

    call VerifyTasks.CompareVCFsVerbosely as CompareFilteredVcfsVerbosely {
      input:
        actual = test_vcfs[idx],
        expected = truth_vcfs[idx],
        actual_index = test_vcf_indexes[idx],
        expected_index = truth_vcf_indexes[idx],
        extra_args = " --ignore-attribute AVERAGE_TREE_SCORE --ignore-attribute CALIBRATION_SENSITIVITY "
      + "--ignore-attribute TREE_SCORE --ignore-attribute SCORE --ignore-filters --ignore-attribute AS_FilterStatus "
      + "--ignore-attribute ExcessHet --ignore-star-attributes --allow-nan-mismatch --ignore-attribute END"
    }
  }

  call MetricsVerification.VerifyMetrics {
    input:
      test_metrics = test_metrics,
      truth_metrics = truth_metrics
  }

  call VerifyTasks.CompareTextFiles as CompareIntervals {
    input:
      test_text_files = test_intervals,
      truth_text_files = truth_intervals
  }

  call CompareFingerprints {
    input:
      test_fingerprint = test_fingerprint,
      truth_fingerprint = truth_fingerprint
  }

  call VerifyNA12878.VerifyNA12878 {
    input:
      vcf_files = flatten([test_vcfs, truth_vcfs]),
      vcf_file_indexes = flatten([test_vcf_indexes, truth_vcf_indexes]),
      vcf_names = ["test","truth"]
  }

  output {
    Array[File]? tables = VerifyNA12878.concordance_tables
  }
  
  meta {
    allowNestedInputs: true
  }
}

task CompareFingerprints {
  input {
    File test_fingerprint
    File truth_fingerprint
  }

  command {
    # We need to strip out the vcf paths because they contain workflow-specific path elements
    cmp <(grep -v '^#' ~{test_fingerprint} | 's/gs:\/\/.*\.vcf\.gz//g') <(grep -v '^#' ~{truth_fingerprint} | 's/gs:\/\/.*\.vcf\.gz//g')
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}