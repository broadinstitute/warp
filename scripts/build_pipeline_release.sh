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

  # Store the current working directory
  CWD=$(pwd)

  outputDir=${OUTPUT_DIR:-${PWD}}/${wdlBasename}
  
  local version=${VERSION}
  if [[ -n "${version}" ]]; then
    version="_${version}"
  fi

  if [[ ! -d ${outputDir} ]]; then
    mkdir -p ${outputDir}
  fi
  # make the outputPrefix use the absolute path of outputDir
  local outputPrefix=$(cd $( dirname $outputDir) && pwd)/$(basename $outputDir)/${wdlBasename}

  local -r outputVersionedPrefix=${outputPrefix}${version}
  
  # Strip the paths out of the root WDL imports
  sed -E '/http/! s/import "(.*)\/(.*\'${WDL_SUFFIX}')"/import "\2"/g' ${rootWdl} > ${outputVersionedPrefix}${WDL_SUFFIX}

  write_options ${rootWdl} ${outputVersionedPrefix}

  versioned_dependencies_zip=${outputVersionedPrefix}${ZIP_SUFFIX}

  cd ${REPO_ROOT}
  write_dependencies_zip ${rootWdl} ${versioned_dependencies_zip}
  cd ${CWD}

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
  local -r rootWdl=${1} versioned_dependencies_zip=${2} working_dir=$(mktemp -d)
  local -r -a dependencies=($(get_dependencies_for_wdl ${rootWdl} | xargs -n1 | sort -u | xargs))

  for file in ${dependencies[@]}; do
    flattened_name=$(basename ${file})
    sed -E '/http/! s/import "(.*)\/(.*\'${WDL_SUFFIX}')"/import "\2"/g' ${file} > ${working_dir}/${flattened_name}
    zip -j ${versioned_dependencies_zip} ${working_dir}/${flattened_name}
  done
  rm -rf ${working_dir}
}

function show_help() {
    echo "Usage: $0 [arguments]"
    echo "This program writes the wdl, dependencies zip and options file for a specified workflow"
    echo ""
    echo "Arguments:"
    echo "  -w             The path to the workflow (.wdl) file"
    echo "  -v             The version of the workflow (used in building the release name)"
    echo "  -o             The directory into which the outputs will be written"
    echo "  -e             The environment (dev, staging, or prod)"
    echo "  -h             print this helpful message"
    echo ""
}

while getopts "hw:v:o:e:" opt; do
    case ${opt} in
	    h)
	      show_help
	      exit 0
	      ;;
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