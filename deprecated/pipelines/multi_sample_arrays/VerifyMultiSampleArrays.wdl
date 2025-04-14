version 1.0

import "../verification/VerifyTasks.wdl" as Tasks

## Copyright Broad Institute, 2018
##
## This WDL script is designed to verify (compare) the outputs of an MultiSampleArrays wdl.
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

workflow VerifyMultiSampleArrays {
  input {
    File truth_vcf
    File test_vcf
    Boolean? done
  }

  call Tasks.CompareVcfsAllowingQualityDifferences {
    input:
      file1 = truth_vcf,
      file2 = test_vcf
  }
  meta {
    allowNestedInputs: true
  }
}
