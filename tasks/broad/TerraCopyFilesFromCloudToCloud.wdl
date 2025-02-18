version 1.0

## Copyright Broad Institute, 2024
##
## This WDL defines tasks used for moving files from place to place on Terra Platform.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

task TerraCopyFilesFromCloudToCloud {
  input {
    Array[String] files_to_copy
    String destination_cloud_path
    Float? contamination
  }

  command {
    set -euo pipefail

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2
    echo ~{default='no_contamination' contamination} > contamination

    if ! grep -q no_contamination contamination; then
      gcloud storage cp -L cp.log contamination ~{destination_cloud_path}.contamination
    fi
    gcloud storage cp ~{sep=' ' files_to_copy} ~{destination_cloud_path}
  }

  output {
    Boolean done = true
  }

  runtime {
    memory: "16 GiB"
    cpu: "1"
    disks: "local-disk 32 HDD"
    docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:499.0.0-slim"
    preemptible: 3
  }
}
