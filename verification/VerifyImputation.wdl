version 1.0

import "../verification/VerifyTasks.wdl" as Tasks
import "../tasks/broad/ImputationTasks.wdl" as ImputationTasks

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
    File haplotype_database
    Boolean split_output_to_single_sample
    String output_callset_name

    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_vcf
    File test_vcf
    File test_vcf_index

    File? input_multi_sample_vcf
    File? input_multi_sample_vcf_index
    Array[File]? input_single_sample_vcfs
    Array[File]? input_single_sample_vcfs_indices

    Array[File]? single_sample_truth_vcf
    Array[File]? single_sample_test_vcf
    Array[File]? single_sample_test_vcf_indices
  }

  scatter (idx in range(length(truth_metrics))) {
    call CompareImputationMetrics {
      input:
        test_metrics = test_metrics[idx],
        truth_metrics = truth_metrics[idx]
    }
  }

  call Tasks.CompareVcfs as CompareOutputVcfs {
    input:
      file1 = truth_vcf,
      file2 = test_vcf
  }

  if (defined(single_sample_truth_vcf)) {
    scatter (idx in range(length(select_first([single_sample_truth_vcf])))) {
      call Tasks.CompareVcfs as CompareSingleSampleOutputVcfs {
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
      haplotypeDatabase = haplotype_database,
      basename = output_callset_name
  }

  if (split_output_to_single_sample) {
    call ImputationTasks.SplitMultiSampleVcf {
      input:
        multiSampleVcf = test_vcf
    }

    call CrosscheckFingerprints as CrosscheckFingerprintsSplit {
      input:
        firstInputs = if (defined(input_multi_sample_vcf)) then select_all([input_multi_sample_vcf]) else select_first([input_single_sample_vcfs]),
        firstInputIndices = if (defined(input_multi_sample_vcf_index)) then select_all([input_multi_sample_vcf_index]) else select_first([input_single_sample_vcfs_indices]),
        secondInputs = SplitMultiSampleVcf.single_sample_vcfs,
        secondInputIndices = SplitMultiSampleVcf.single_sample_vcf_indices,
        haplotypeDatabase = haplotype_database,
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
    File haplotypeDatabase
    String basename
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
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

    gatk CrosscheckFingerprints -I ~{sep=" -I " firstInputs} -SI ~{sep=" -SI " secondInputs} -H ~{haplotypeDatabase} -O ~{basename}.crosscheck
  >>>

  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " HDD"
    memory: "8 GiB"
  }

  output {
    File crosscheck = "~{basename}.crosscheck"
  }
}