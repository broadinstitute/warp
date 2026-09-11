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
    # popped posteriors vcf, one per chromosome
    Array[File] truth_imputed_vcf
    Array[File] test_imputed_vcf

    Boolean? done
  }


  scatter (idx in range(length(truth_imputed_vcf))) {
    call Tasks.CompareVcfs as CompareImputedVcfs {
      input:
        file1 = truth_imputed_vcf[idx],
        file2 = test_imputed_vcf[idx],
        patternForLinesToExcludeFromComparison = "##"
    }
  }

  output {
  }
  meta {
    allowNestedInputs: true
  }
}
