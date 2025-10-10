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

workflow VerifyArrayImputationQC {
  input {
    File truth_outputs
    File test_outputs

    Boolean? done
  }

  call CompareFiles {
    input:
      truth_outputs = truth_outputs,
      test_outputs = test_outputs
  }

  output {
  }
}

task CompareFiles {
  input {
    File truth_outputs
    File test_outputs
  }
  command <<<
    set -eo pipefail
    diff "~{test_outputs}" "~{truth_outputs}"

    if [ $? -ne 0 ];
    then
      echo "Error: ${test_outputs} and ${truth_outputs} differ"
    fi
  >>>

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}
