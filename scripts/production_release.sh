#!/usr/bin/env bash

set -e

declare -r SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
source ${SCRIPT_DIR}/common.sh

declare -r PREVIOUS_COMMIT="HEAD~1"

${SCRIPT_DIR}/validate_release.sh -g ${PREVIOUS_COMMIT}

pipelines_to_release=($(get_modified_pipelines ${PREVIOUS_COMMIT}))

for pipeline in ${pipelines_to_release[@]}; do
  ${SCRIPT_DIR}/release_pipeline_to_github.sh -p ${pipeline} -e prod
done
