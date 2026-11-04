version 1.0

## Copyright Broad Institute, 2010
##
## This WDL defines tasks used for moving files from place to place on Google Cloud Platform.
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

task CopyFilesFromCloudToCloud {
  input {
    Array[String] files_to_copy
    String destination_cloud_path
    String vault_token_path
    String google_account_vault_path
    Float? contamination
    String base_file_name = "base_file"
  }

  command {
    set -euo pipefail
    export PATH=$PATH:/usr/local/google-cloud-sdk/bin

    export VAULT_ADDR=https://clotho.broadinstitute.org:8200
    export VAULT_TOKEN=$(gsutil cat ~{vault_token_path})

    vault read -format=json ~{google_account_vault_path} | jq .data > picard-account.pem
    /usr/local/google-cloud-sdk/bin/gcloud auth activate-service-account --key-file=picard-account.pem


    echo ~{default='no_contamination' contamination} > contamination

    # We use gsutil here to copy files to a potentially different cloud location - in parallel.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf tmp_dest
    mkdir tmp_dest
      RETRY_LIMIT=5

    count=0
    until cat ~{write_lines(files_to_copy)} | /usr/local/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I ~{destination_cloud_path}; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if ! grep -q no_contamination contamination; then
      /usr/local/google-cloud-sdk/bin/gsutil -m cp -L cp.log contamination ~{destination_cloud_path}~{base_file_name}.contamination
    fi
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the files to the cloud destination' && exit 1
    fi
  }

  output {
    Boolean done = true
  }

  # The 'noAddress' runtime parameter is set to false here because
  # Vault needs to talk to the Broad Vault server to get auth information.
  # In the future, we should store the extracted data in a GCS bucket so that
  # We don't have to use an external IP.
  runtime {
    memory: "2 GiB"
    cpu: "1"
    disks: "local-disk 20 HDD"
    docker: "us.gcr.io/broad-gotc-prod/dsde-toolbox:stable_04-18-2022"
    preemptible: 3
    noAddress: false
  }
}
