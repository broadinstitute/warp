#!/usr/bin/env bash
set -e -o pipefail

declare -r SCRIPT_DIR=$(cd $(dirname $0) && pwd)
declare -r REPO_ROOT=$(cd $(dirname ${SCRIPT_DIR}) && pwd)

declare -r -a ZIP_DIRS=(pipelines beta-pipelines structs tasks verification tests/skylab)
declare -r ZIP_PREFIX=workflow_dependencies


function deploy_dependencies() {
  local -r prefix="$1" version="$2" target_dir="$3" preserve_dir_structure="$4"

  local -r versioned_dependencies_zip=${prefix}.${version}.zip
  local -r dependencies_zip_link=${prefix}.zip

  if [[ ${preserve_dir_structure} == 'preserve_dir_structure' ]]; then
    cd ${REPO_ROOT}
    zip -r ${target_dir}/${versioned_dependencies_zip} ${ZIP_DIRS[@]}
  else
    local -r working_dir=$(mktemp -d)
    cd ${REPO_ROOT}
    for file in $(find ${ZIP_DIRS[@]} -type f -name '*.wdl'); do
      flattened_name=$(basename ${file})
      sed -E 's/import "(.*)\/(.*\.wdl)"/import "\2"/g' ${file} > ${working_dir}/${flattened_name}
      zip -u -j -v ${target_dir}/${versioned_dependencies_zip} ${working_dir}/${flattened_name}
    done
    rm -rf ${working_dir}
  fi
  cd ${target_dir}
  ln -sfv ${versioned_dependencies_zip} ${dependencies_zip_link}
  cd ${REPO_ROOT}
}

function deploy_options() {
  local -r prefix="$1" version="$2" wdl_dir="$3" target_dir="$4" env="$5"

  local -r versioned_options=${prefix}.${version}.options.json

  # Some workflows have per-environment options, and others share options across all environments.
  local base_options
  if [ -f ${wdl_dir}/${prefix}.options.json ]; then
    base_options=${wdl_dir}/${prefix}.options.json
  elif [ -f ${wdl_dir}/${prefix}.${env}.options.json ]; then
    base_options=${wdl_dir}/${prefix}.${env}.options.json
  elif [[ "${wdl_dir}" =~ .*"skylab".* ]] && [[ -f ${REPO_ROOT}/tests/skylab/test.options.json ]]; then
    base_options=${REPO_ROOT}/tests/skylab/test.options.json
  else
    echo >&2 Error: Options JSON not found at either ${prefix}.options.json or ${prefix}.${env}.options.json
    exit 1
  fi

  cp ${base_options} ${target_dir}/${versioned_options}
  cd ${target_dir}
  ln -sfv ${versioned_options} ${prefix}.options.json
  cd ${REPO_ROOT}
}

function deploy_wdl() {
  local -r prefix="$1" version="$2" wdl_dir="$3" target_dir="$4" preserve_dir_structure="$5"

  local -r base_wdl=${prefix}.wdl
  local -r versioned_wdl=${prefix}.${version}.wdl

  cp ${wdl_dir}/${base_wdl} ${target_dir}/${versioned_wdl}
  local sed_command=''
  if [[ ${preserve_dir_structure} == 'preserve_dir_structure' ]]; then
    sed_command='s/import "(\.\.\/)*(.*)"/import "\2"/g'
  else
    sed_command='s/import "(.*)\/(.*\.wdl)"/import "\2"/g'
  fi
    sed -i.unedited -E "${sed_command}" ${target_dir}/${versioned_wdl} && rm ${target_dir}/${versioned_wdl}.unedited
  cd ${target_dir}
  ln -sfv ${versioned_wdl} ${base_wdl}
  cd ${REPO_ROOT}
}

function main() {
    local -r wdl="$1" cloud_workflows_target="$2" pipeline_hash="$3" env="$4" preserve_dir_structure="$5"
    local -r pwd=$(pwd)
    trap "cd ${pwd}" ERR EXIT HUP INT TERM

    if [ ! -f ${wdl} ]; then
      echo >&2 Error: ${wdl} does not exist!
      exit 1
    fi

    local -r wdl_basename=${wdl##*/}
    local -r wdl_prefix=${wdl_basename%.*}
    local -r wdl_dir=$(cd $(dirname ${wdl}) && pwd)

    local -r deps_dir=${cloud_workflows_target}/dependencies
    local -r pipeline_dir=${cloud_workflows_target}/${wdl_prefix}

    mkdir -p ${deps_dir}
    mkdir -p ${pipeline_dir}

    deploy_dependencies ${ZIP_PREFIX} ${pipeline_hash} ${deps_dir} ${preserve_dir_structure}
    deploy_options ${wdl_prefix} ${pipeline_hash} ${wdl_dir} ${pipeline_dir} ${env}
    deploy_wdl ${wdl_prefix} ${pipeline_hash} ${wdl_dir} ${pipeline_dir} ${preserve_dir_structure}
}

main "${@}"
