#!/usr/bin/env bash

declare -r SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
source ${SCRIPT_DIR}/common.sh

declare ROOT_WDL=""
declare VERSION=""
declare OUTPUT_DIR=""
declare ENV=""

declare -r ZIP_SUFFIX=".zip" WDL_SUFFIX=".wdl" OPTIONS_SUFFIX=".options.json"

function make_release() {
  local -r rootWdl=${ROOT_WDL}
  local -r wdlBasename=$(basename ${rootWdl} ${WDL_SUFFIX})

  outputDir=${OUTPUT_DIR:-${PWD}}/${wdlBasename}
  
  local version=${VERSION}
  if [[ -n "${version}" ]]; then
    version="_${version}"
  fi

  if [[ ! -d ${outputDir} ]]; then
    mkdir -p ${outputDir}
  fi

  local outputPrefix=${outputDir}/${wdlBasename}

  local -r outputVersionedPrefix=${outputPrefix}${version}
  
  # Strip the paths out of the root WDL imports
  sed -E 's/import "(.*)\/(.*\'${WDL_SUFFIX}')"/import "\2"/g' ${rootWdl} > ${outputVersionedPrefix}${WDL_SUFFIX}

  write_options ${rootWdl} ${outputVersionedPrefix}
  write_dependencies_zip ${rootWdl} ${outputVersionedPrefix}

}

function write_options() {
  local -r rootWdl=${1} versioned_options=${2}${OPTIONS_SUFFIX}
  local -r rootWdlDir=$(dirname ${rootWdl})
  local -r rootWdlName=$(basename ${rootWdl} ${WDL_SUFFIX})
  
  local env=${ENV}
  if [[ -n "${env}" ]]; then
    env=".${env}"
  fi

  # Some workflows have per-environment options, and others share options across all environments.
  local base_options
  if [[ -f ${rootWdlDir}/${rootWdlName}${OPTIONS_SUFFIX} ]]; then
    base_options=${rootWdlDir}/${rootWdlName}${OPTIONS_SUFFIX}
  elif [[ -f ${rootWdlDir}/${rootWdlName}${env}${OPTIONS_SUFFIX} ]]; then
    base_options=${rootWdlDir}/${rootWdlName}${env}${OPTIONS_SUFFIX}
  else
    echo "Writing an empty options file"
    tmpOptions=$(mktemp)
    echo '{}' > ${tmpOptions}
    base_options=${tmpOptions}
  fi

  cp ${base_options} ${versioned_options}
}

function write_dependencies_zip() {
  local -r rootWdl=${1} versioned_dependencies_zip=${2}${ZIP_SUFFIX} working_dir=$(mktemp -d)
  local -r -a dependencies=($(get_dependencies_for_wdl ${rootWdl} | xargs -n1 | sort -u | xargs))

  # cd ${DSDE_PIPELINES_ROOT}
  for file in ${dependencies[@]}; do
    flattened_name=$(basename ${file})
    sed -E 's/import "(.*)\/(.*\'${WDL_SUFFIX}')"/import "\2"/g' ${file} > ${working_dir}/${flattened_name}
    zip -j ${versioned_dependencies_zip} ${working_dir}/${flattened_name}
  done
  rm -rf ${working_dir}
}

while getopts "w:v:o:e:" opt; do
    case ${opt} in 
      w)
        if [[ ! -f ${OPTARG} ]]; then
          echo >&2 Error: ${OPTARG} does not exist!
          exit 1
        fi
        ROOT_WDL=${OPTARG}
        ;;
      v)
        VERSION=${OPTARG}
        ;;
      o)
        OUTPUT_DIR=${OPTARG}
        ;;
      e)
        ENV=${OPTARG}
        ;;
      *)
        >&2 echo "invalid option"
        exit 1
        ;; 
    esac 
  done
  shift $((OPTIND-1)) 

make_release