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
      echo CloudWorkflow -p AnnotationFiltration -t ${test} --env ${env};;
    Arrays)
      echo CloudWorkflow -p Arrays -a Single ${common_args};;
    MultiSampleArrays)
      echo CloudWorkflow -p Arrays -a Multi ${common_args};;
    BroadInternalRNAWithUMIs)
      echo CloudWorkflow -p BroadInternalRNAWithUMIs ${common_args};;
    BroadInternalUltimaGenomics)
      echo CloudWorkflow -p BroadInternalUltimaGenomics ${common_args};;
    CheckFingerprint)
      echo CloudWorkflow -p CheckFingerprint ${common_args};;
    ExomeGermlineSingleSample)
      echo CloudWorkflow -p GermlineSingleSample -d Exome ${common_args};;
    ExomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p Reprocessing -d Exome --env ${env} -t Plumbing -b ${truth} ${uncached}
      else
        continue
      fi;;
    JointGenotyping)
      echo CloudWorkflow -p JointGenotyping -d Exome ${common_args} --papi-version PAPIv2;
      echo CloudWorkflow -p JointGenotyping -d WGS --env ${env} -t Plumbing -b ${truth} ${uncached} --papi-version PAPIv2;;
    IlluminaGenotypingArray)
      echo CloudWorkflow -p IlluminaGenotypingArray ${common_args};;
    Imputation)
      echo CloudWorkflow -p Imputation ${common_args};;
    ExternalExomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p ExternalReprocessing -d Exome --env ${env} -t Plumbing -b ${truth} ${uncached}
      else
        continue
      fi;;
    ExternalWholeGenomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p ExternalReprocessing -d WGS --env ${env} -t Plumbing -b ${truth} ${uncached}
      else
        continue
      fi;;
    WholeGenomeGermlineSingleSample)
      echo CloudWorkflow -p GermlineSingleSample -d WGS ${common_args};;
    WholeGenomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p Reprocessing -d WGS --env ${env} -t Plumbing -b ${truth} ${uncached}
      else
        continue
      fi;;
    ValidateChip)
      echo CloudWorkflow -p ValidateChip ${common_args};;
    ReblockGVCF)
      echo CloudWorkflow -p ReblockGvcf -d Exome ${common_args};
      echo CloudWorkflow -p ReblockGvcf -d WGS ${common_args};;
    RNAWithUMIsPipeline)
      echo CloudWorkflow -p RNAWithUMIs ${common_args};;
    TargetedSomaticSingleSample)
      echo CloudWorkflow -p SomaticSingleSample -d Targeted ${common_args};;
    CramToUnmappedBams)
      echo CloudWorkflow -p CramToUnmappedBams ${common_args};;
    JointGenotypingByChromosomePartOne)
      continue;;
    JointGenotypingByChromosomePartTwo)
      continue;;
    UltimaGenomicsGermlineSingleSample)
      echo CloudWorkflow -p UltimaGenomicsGermlineSingleSample ${common_args};;
    UltimaGenomicsJointGenotyping)
      echo CloudWorkflow -p UltimaGenomicsJointGenotyping ${common_args};;
    GDCWholeGenomeSomaticSingleSample)
      echo CloudWorkflow -p GDCWholeGenomeSomaticSingleSample -d WGS ${common_args};;
    VariantCalling)
      echo CloudWorkflow -p VariantCalling -d Exome -t Plumbing --env ${env} -b ${truth} ${uncached};
      echo CloudWorkflow -p VariantCalling -d WGS -t Plumbing --env ${env} -b ${truth} ${uncached};;
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