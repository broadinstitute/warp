version 1.0

import "../verification/VerifyTasks.wdl" as Tasks

## Copyright Broad Institute, 2018
##
## This WDL script is designed to verify (compare) the outputs of Glimpse2 SV Imputation wdl.
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

workflow VerifyGlimpse2SVImputation {
  input {
    # bubble posteriors vcf, one per chromosome
    Array[File] truth_bubble_posteriors_vcf
    Array[File] test_bubble_posteriors_vcf

    # popped posteriors vcf, one per chromosome
    Array[File] truth_popped_posteriors_vcf
    Array[File] test_popped_posteriors_vcf

    Boolean? done
  }

  scatter (idx in range(length(truth_bubble_posteriors_vcf))) {
    call Tasks.CompareVcfs as CompareBubblePosteriorsVcfs {
      input:
        file1 = truth_bubble_posteriors_vcf[idx],
        file2 = test_bubble_posteriors_vcf[idx]
    }
  }

  scatter (idx in range(length(truth_popped_posteriors_vcf))) {
    call Tasks.CompareVcfs as ComparePoppedPosteriorsVcfs {
      input:
        file1 = truth_popped_posteriors_vcf[idx],
        file2 = test_popped_posteriors_vcf[idx]
    }
  }

  output {
  }
  meta {
    allowNestedInputs: true
  }
}
