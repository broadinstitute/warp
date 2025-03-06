version 1.0

import "../verification/VerifyIlluminaGenotypingArray.wdl" as VerifyIlluminaGenotypingArray
import "../verification/VerifyTasks.wdl" as VerifyTasks

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

workflow VerifyArrays {
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

    File truth_params_file
    File test_params_file

    Boolean? done
  }

  call VerifyIlluminaGenotypingArray.VerifyIlluminaGenotypingArray {
    input:
      truth_metrics = truth_metrics,
      test_metrics = test_metrics,
      test_gtc = test_gtc,
      truth_gtc = truth_gtc,
      bead_pool_manifest_file = bead_pool_manifest_file,
      truth_vcf = truth_vcf,
      test_vcf = test_vcf,
      truth_fp_vcf = truth_fp_vcf,
      test_fp_vcf = test_fp_vcf,
      truth_red_idat_md5 = truth_red_idat_md5,
      test_red_idat_md5 = test_red_idat_md5,
      truth_green_idat_md5 = truth_green_idat_md5,
      test_green_idat_md5 = test_green_idat_md5
  }

  call VerifyTasks.CompareTextFiles as CompareParamsFiles {
    input:
      test_text_files = [test_params_file],
      truth_text_files = [truth_params_file]
  }

  output {
  }
  meta {
    allowNestedInputs: true
  }
}

