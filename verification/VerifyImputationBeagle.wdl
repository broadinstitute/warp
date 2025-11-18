version 1.0

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

workflow VerifyImputationBeagle {
  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    # imputed variant multi sample vcf
    File multi_sample_truth_vcf
    File multi_sample_test_vcf

    # imputed hom ref sites only vcf
    File sites_only_truth_vcf
    File sites_only_test_vcf

    Boolean? done
  }

  String bcftools_docker_tag = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"

  scatter (idx in range(length(truth_metrics))) {
    call CompareImputationMetrics {
      input:
        test_metrics = test_metrics[idx],
        truth_metrics = truth_metrics[idx]
    }
  }

  call Tasks.CompareVcfs as CompareOutputMultiSampleVcfs {
    input:
      file1 = multi_sample_truth_vcf,
      file2 = multi_sample_test_vcf,
      patternForLinesToExcludeFromComparison = "##" # ignore headers
  }

  call Tasks.CompareVcfs as CompareOutputSitesOnlyVcfs {
    input:
      file1 = sites_only_truth_vcf,
      file2 = sites_only_test_vcf,
      patternForLinesToExcludeFromComparison = "##" # ignore headers
  }

  output {
  }
  meta {
    allowNestedInputs: true
  }
}

task CompareImputationMetrics {
  input {
    File test_metrics
    File truth_metrics
  }
  command <<<
    set -eo pipefail
    diff "~{test_metrics}" "~{truth_metrics}"

    if [ $? -ne 0 ];
    then
      echo "Error: ${test_metrics} and ${truth_metrics}  differ"
    fi
  >>>

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}
