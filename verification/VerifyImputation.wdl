version 1.0

import "../verification/VerifyTasks.wdl" as Tasks
import "../tasks/wdl/ImputationTasks.wdl" as ImputationTasks

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

workflow VerifyImputation {
  input {
    Boolean split_output_to_single_sample
    String output_callset_name

    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_vcf
    File test_vcf
    File test_vcf_index
    File truth_vcf_index

    File? input_multi_sample_vcf
    File? input_multi_sample_vcf_index
    Array[File]? input_single_sample_vcfs
    Array[File]? input_single_sample_vcfs_indices

    Array[File]? single_sample_truth_vcf
    Array[File]? single_sample_test_vcf

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

  call Tasks.CompareVcfsAllowingQualityDifferences as CompareOutputVcfs {
    input:
      file1 = truth_vcf,
      file2 = test_vcf
  }

  if (defined(single_sample_truth_vcf)) {
    scatter (idx in range(length(select_first([single_sample_truth_vcf])))) {
      call Tasks.CompareVcfsAllowingQualityDifferences as CompareSingleSampleOutputVcfs {
        input:
          file1 = select_first([single_sample_test_vcf])[idx],
          file2 = select_first([single_sample_truth_vcf])[idx]
      }
    }
  }

  call CrosscheckFingerprints {
    input:
      firstInputs = if (defined(input_multi_sample_vcf)) then select_all([input_multi_sample_vcf]) else select_first([input_single_sample_vcfs]),
      firstInputIndices = if (defined(input_multi_sample_vcf_index)) then select_all([input_multi_sample_vcf_index]) else select_first([input_single_sample_vcfs_indices]),
      secondInputs = [test_vcf],
      secondInputIndices = [test_vcf_index],
      basename = output_callset_name
  }

  if (split_output_to_single_sample) {
    call ImputationTasks.CountSamples {
      input:
        vcf = test_vcf
    }
    call ImputationTasks.SplitMultiSampleVcf {
      input:
        multiSampleVcf = test_vcf,
        bcftools_docker = bcftools_docker_tag,
        nSamples = CountSamples.nSamples
    }

    call CrosscheckFingerprints as CrosscheckFingerprintsSplit {
      input:
        firstInputs = if (defined(input_multi_sample_vcf)) then select_all([input_multi_sample_vcf]) else select_first([input_single_sample_vcfs]),
        firstInputIndices = if (defined(input_multi_sample_vcf_index)) then select_all([input_multi_sample_vcf_index]) else select_first([input_single_sample_vcfs_indices]),
        secondInputs = SplitMultiSampleVcf.single_sample_vcfs,
        secondInputIndices = SplitMultiSampleVcf.single_sample_vcf_indices,
        basename = output_callset_name
    }
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

task CrosscheckFingerprints {
  input {
    Array[File] firstInputs
    Array[File] secondInputs
    Array[File] firstInputIndices
    Array[File] secondInputIndices
    String basename
    File haplotypeDatabase = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.haplotype_database.txt"
    String picard_docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
  }

  Int disk_size = ceil(1.2*(size(firstInputs, "GiB") + size(secondInputs, "GiB") + size(haplotypeDatabase, "GiB"))) + 100

  command <<<
    # add links to ensure correctly located indices
    array_vcfs=( ~{sep=" " firstInputs} )
    array_indices=( ~{sep=" " firstInputIndices} )
    for i in ${!array_vcfs[@]}; do
      ln -s ${array_indices[i]} $(dirname ${array_vcfs[i]})
    done

    array_vcfs2=( ~{sep=" " secondInputs} )
    array_indices2=( ~{sep=" " secondInputIndices} )
    for i in ${!array_vcfs2[@]}; do
      ln -s ${array_indices2[i]} $(dirname ${array_vcfs2[i]})
    done

    java -Xms3500m -Xms7500m -jar /usr/picard/picard.jar \
      CrosscheckFingerprints \
      -I ~{sep=" -I " firstInputs} \
      -SI ~{sep=" -SI " secondInputs} \
      -H ~{haplotypeDatabase} \
      -O ~{basename}.crosscheck
  >>>

  runtime {
    docker: picard_docker
    disks: "local-disk " + disk_size + " HDD"
    memory: "8000 MiB"
  }

  output {
    File crosscheck = "~{basename}.crosscheck"
  }
}