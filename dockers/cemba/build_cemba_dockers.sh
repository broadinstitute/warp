#!/usr/bin/env bash
set -euo pipefail

# relative location of additonal scripts
declare -r SCRIPT_DIR="/scripts"

# versions of tools used
declare -r BOWTIE2_VERSION=2.3.4.3
declare -r BISMARK_VERSION=0.21.0
declare -r CUTADAPT_VERSION=1.18
declare -r SAMTOOLS_VERSION=1.9
declare -r PICARD_VERSION=2.18.23
declare -r BISULFITE_REF_VERSION=1.0
declare -r MONITORING_VERSION=1.0

function build_docker () {
  # using the tool/dir name
  local -r image_tag=$1
  # get the version of the tool/docker
  local -r image_version=$(echo $2 | cut -d"=" -f2)
  # get the relative path to the docker
  local -r dockerfile_path="docker/${image_tag}/Dockerfile"
  # get the images name
  local -r image_name="quay.io/broadinstitute/${image_tag}:${image_version}"
  # get additional arguments with the "--build-arg" option
  local -a build_args=()
  for arg in $@; do
    build_args=(--build-arg ${arg} ${build_args[@]})
  done

  # build and push the docker
  docker build ${build_args[@]} -t ${image_name} -f ${dockerfile_path} .
  docker push ${image_name}
}

# build dockers (dir/tool name must be 1st arg & version name must be 2nd arg)
build_docker samtools SAMTOOLS_VERSION=${SAMTOOLS_VERSION}
build_docker bowtie2 BOWTIE2_VERSION=${BOWTIE2_VERSION}
build_docker bismark BISMARK_VERSION=${BISMARK_VERSION}
build_docker bisulfite-references BISULFITE_REF_VERSION=${BISULFITE_REF_VERSION} SCRIPT_DIR=${SCRIPT_DIR}
build_docker cutadapt CUTADAPT_VERSION=${CUTADAPT_VERSION}
build_docker picard PICARD_VERSION=${PICARD_VERSION}
build_docker monitoring MONITORING_VERSION=${MONITORING_VERSION} SCRIPT_DIR=${SCRIPT_DIR}

