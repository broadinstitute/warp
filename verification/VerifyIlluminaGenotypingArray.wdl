version 1.0

import "../verification/VerifyMetrics.wdl" as MetricsVerification
import "../verification/VerifyTasks.wdl" as Tasks

## Copyright Broad Institute, 2018
##
## This WDL script is designed to verify (compare) the outputs of an ArrayWf wdl.
##
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow VerifyIlluminaGenotypingArray {
  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    File test_gtc
    File truth_gtc
    File bead_pool_manifest_file

    File truth_vcf
    File test_vcf

    File truth_fp_vcf
    File test_fp_vcf

    File truth_red_idat_md5
    File test_red_idat_md5

    File truth_green_idat_md5
    File test_green_idat_md5

    Boolean? done
  }

  call MetricsVerification.CompareTwoNumbers {
    input:
      num1 = length(test_metrics),
      num2 = length(truth_metrics),
      error_msg = "Different number of metric files"
  }

  String bafm_ext = ".bafregress_metrics"
  String avcdm_ext = "arrays_variant_calling_detail_metrics"

  scatter (idx in range(length(truth_metrics))) {
    String metrics_basename = basename(truth_metrics[idx])
    Boolean is_bafm_file = basename(metrics_basename, bafm_ext) != metrics_basename
    if (!is_bafm_file) {
      Boolean is_avcdm_file = basename(metrics_basename, avcdm_ext) != metrics_basename
      Array[String] metrics_to_ignore = if (is_avcdm_file) then ["AUTOCALL_DATE", "ANALYSIS_VERSION", "PIPELINE_VERSION"] else []
      call MetricsVerification.CompareMetricFiles {
        input:
          dependency_input = CompareTwoNumbers.output_file,
          file1 = test_metrics[idx],
          file2 = truth_metrics[idx],
          output_file = "metric_~{idx}.txt",
          metrics_to_ignore = metrics_to_ignore
      }
    }
    if (is_bafm_file) {
      # TODO - is this still necessary??
      # BAF metrics files are defined in picard-private.  The tool used by CompareMetricFiles (above) cannot
      # instantiate them and so fails.  This is a hopefully short-term fix for that problem - just compare the metrics files as text files.
      call CompareMetricFilesAsText {
        input:
          file1 = test_metrics[idx],
          file2 = truth_metrics[idx]
      }
    }
  }

  call Tasks.CompareVcfsAllowingQualityDifferences as CompareOutputVcfs {
    input:
      file1 = truth_vcf,
      file2 = test_vcf
  }

  call Tasks.CompareVcfsAllowingQualityDifferences as CompareOutputFingerprintVcfs {
    input:
       file1 = truth_fp_vcf,
       file2 = test_fp_vcf
  }

  call Tasks.CompareGtcs {
    input:
      file1 = test_gtc,
      file2 = truth_gtc,
      bead_pool_manifest_file = bead_pool_manifest_file
  }

  call CompareFiles as CompareGreenIdatMdf5Sum {
    input:
      file1 = test_green_idat_md5,
      file2 = truth_green_idat_md5
  }

  call CompareFiles as CompareRedIdatMdf5Sum {
    input:
      file1 = test_red_idat_md5,
      file2 = truth_red_idat_md5
  }
  output {
  }
  meta {
    allowNestedInputs: true
  }
}

task CompareGtcs {
  input {
    File file1
    File file2
    File illumina_normalization_manifest
  }

  command {
    java -Xmx3g -Dpicard.useLegacyParser=false  -jar /usr/picard/picard.jar \
      CompareGtcFiles \
      --INPUT ~{file1} \
      --INPUT ~{file2} \
      --ILLUMINA_NORMALIZATION_MANIFEST ~{illumina_normalization_manifest}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
  }
}

task CompareFiles {

  input {
    File file1
    File file2
  }

  command {
    if ! diff ~{file1} ~{file2}; then
       echo "Error: Files ~{file1} and ~{file2} differ" >&2
       exit 1
    fi
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}

task CompareMetricFilesAsText {

  input {
    File file1
    File file2
  }

  command {
    # Compares two text files, ignoring lines that begin with a '#'
    diff <(grep -v '# ' ~{file1}) <(grep -v '# ' ~{file2})
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}
