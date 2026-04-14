version 1.0

## Copyright Broad Institute, 2021
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3)
## (see LICENSE in https://github.com/openwdl/wdl).

task pytorch_cuda_test {

  input {
    # Docker image for cellbender remove-background version
    String docker_image

    # Hardware-related inputs
    String? hardware_zones = "us-central1-a us-central1-c"

    Int? hardware_disk_size_GB = 50
    Int? hardware_boot_disk_size_GB = 20
    Int? hardware_cpu_count = 4
    Int? hardware_memory_GB = 15
    String? hardware_gpu_type = "nvidia-tesla-t4"
    String? nvidiaDriverVersion = "418.87.00"  # cromwell default
  }

  command {
    set -e

    echo "which python"
    echo $(which python)
    echo "Starting python code ..."

    python <<CODE

    import torch
    import sys

    print("Imported pytorch")

    print("torch.cuda.is_available()")
    print(torch.cuda.is_available())

    if not torch.cuda.is_available():
        sys.exit(1)

    CODE
  }

  output {
    File log = stdout()
  }

  runtime {
    docker: "${docker_image}"
    bootDiskSizeGb: hardware_boot_disk_size_GB
    disks: "local-disk ${hardware_disk_size_GB} HDD"
    memory: "${hardware_memory_GB}G"
    cpu: hardware_cpu_count
    zones: "${hardware_zones}"
    gpuCount: 1
    gpuType: "${hardware_gpu_type}"
    nvidiaDriverVersion: "${nvidiaDriverVersion}"
    maxRetries: 0
  }

}

workflow test_cuda {

  call pytorch_cuda_test

  output {
    File log = pytorch_cuda_test.log
  }

}
