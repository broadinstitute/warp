#!/usr/bin/env bash

declare SCRIPT_DIR
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
source ${SCRIPT_DIR}/common.sh

declare ALL_PIPELINES=($(get_versioned_pipelines))

function pipeline_to_args() {
  local -r pipeline=${1}
  local -r env=${2}
  local -r test=${3}
  local -r truth=${4}
  local -r uncached=${5}

  local -r common_args="--env ${env} -t ${test} -b ${truth} ${uncached}"

  case ${pipeline} in
    AnnotationFiltration)
      echo AnnotationFiltration -t ${test} --env ${env};;
    Arrays)
      echo Arrays -a Single ${common_args};;
    MultiSampleArrays)
      echo Arrays -a Multi ${common_args};;
    ExomeGermlineSingleSample)
      echo GermlineSingleSample -d Exome ${common_args};;
    ExomeReprocessing)
      echo Reprocessing -d Exome ${common_args};;
    JointGenotyping)
      echo JointGenotyping -d Exome ${common_args} --papi-version PAPIv2;
      echo JointGenotyping -d WGS --env ${env} -t Plumbing -b ${truth} ${uncached} --papi-version PAPIv2;;
    IlluminaGenotypingArray)
      echo IlluminaGenotypingArray ${common_args};;
    ExternalExomeReprocessing)
      echo ExternalReprocessing -d Exome --env ${env} -t Plumbing -b ${truth} ${uncached};;
    ExternalWholeGenomeReprocessing)
      echo ExternalReprocessing -d WGS --env ${env} -t Plumbing -b ${truth} ${uncached};;
    WholeGenomeGermlineSingleSample)
      echo GermlineSingleSample -d WGS ${common_args};;
    WholeGenomeReprocessing)
      echo Reprocessing -d WGS ${common_args};;
    ValidateChip)
      echo ValidateChip ${common_args};;
    ReblockGVCF)
      echo ReblockGvcf -d Exome ${common_args};
      echo ReblockGvcf -d WGS ${common_args};;
    TargetedSomaticSingleSample)
      echo SomaticSingleSample -d Targeted ${common_args};;
    JointGenotypingByChromosomePartOne)
      continue;;
    JointGenotypingByChromosomePartTwo)
      continue;;
  esac
}

function main() {
  local -r gittish=${1}
  local -r test_all=${2}
  local -r env=${3}
  local -r test=${4}
  local -r truth=${5}
  local -r uncached=${6}

  local -a changed_pipeline_paths=()
  local -a args=()

  if ${test_all}; then
    changed_pipeline_paths=(${ALL_PIPELINES[@]})
  else
    changed_pipeline_paths=("$(get_pipelines_to_test ${gittish})")
  fi

  for changed_pipeline_path in ${changed_pipeline_paths[*]}; do
    pipeline=$(basename ${changed_pipeline_path} .wdl)
    arg="$(pipeline_to_args ${pipeline} ${env} ${test} ${truth} ${uncached})"
    if [[ -n ${arg} ]]; then
      args+=("${arg}")
    fi
  done

  for arg in "${args[@]}"; do
    echo "${arg}"
  done
}

main ${1} ${2} ${3} ${4} ${5} ${6}
