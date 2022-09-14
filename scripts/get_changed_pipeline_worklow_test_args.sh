#!/usr/bin/env bash

declare SCRIPT_DIR
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
source ${SCRIPT_DIR}/common.sh

declare ALL_PIPELINES=($(get_versioned_pipelines))

function pipeline_to_args() {
  local -r pipeline=${1}
  local -r test=${2}

  local -r common_args="-t ${test}"

  case ${pipeline} in
    AnnotationFiltration)
      continue;;
    Arrays)
      echo CloudWorkflow -p Arrays ${common_args};;
    BroadInternalRNAWithUMIs)
      echo CloudWorkflow -p BroadInternalRNAWithUMIs ${common_args};;
    BroadInternalUltimaGenomics)
      echo CloudWorkflow -p BroadInternalUltimaGenomics ${common_args};;
    # CEMBA)
    #   echo CloudWorkflow -p CEMBA ${common_args};;
    CheckFingerprint)
      echo CloudWorkflow -p CheckFingerprint ${common_args};;
    CramToUnmappedBams)
      echo CloudWorkflow -p CramToUnmappedBams ${common_args};;
    ExternalExomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p ExternalExomeReprocessing -t Plumbing
      else
        continue
      fi;;
    ExomeGermlineSingleSample)
      echo CloudWorkflow -p ExomeGermlineSingleSample ${common_args};;
    ExomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p ExomeReprocessing -t Plumbing
      else
        continue
      fi;;
    ExternalWholeGenomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p ExternalWholeGenomeReprocessing -t Plumbing
      else
        continue
      fi;;
    GDCWholeGenomeSomaticSingleSample)
      echo CloudWorkflow -p GDCWholeGenomeSomaticSingleSample ${common_args};;
    IlluminaGenotypingArray)
      echo CloudWorkflow -p IlluminaGenotypingArray ${common_args};;
    Imputation)
      echo CloudWorkflow -p Imputation ${common_args};;
    JointGenotyping)
      echo CloudWorkflow -p JointGenotyping ${common_args};;
    JointGenotypingByChromosomePartOne)
      continue;;
    JointGenotypingByChromosomePartTwo)
      continue;;
    MultiSampleArrays)
      echo CloudWorkflow -p MultiSampleArrays ${common_args};;
    MultiSampleSmartSeq2)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p MultiSampleSmartSeq2 -t Plumbing
      else
        echo CloudWorkflow -p MultiSampleSmartSeq2 ${common_args}
      fi;;
    MultiSampleSmartSeq2SingleNucleus)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p MultiSampleSmartSeq2SingleNucleus -t Plumbing
      else
        echo CloudWorkflow -p MultiSampleSmartSeq2SingleNucleus ${common_args}
      fi;;
    Optimus)
      echo CloudWorkflow -p Optimus ${common_args};;
    ReblockGVCF)
      echo CloudWorkflow -p ReblockGvcf  ${common_args};;
    RNAWithUMIsPipeline)
      echo CloudWorkflow -p RNAWithUMIsPipeline ${common_args};;
    scATAC)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p scATAC -t Plumbing
      else
        echo CloudWorkflow -p scATAC ${common_args}
      fi;;
    SmartSeq2SingleSample)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p SmartSeq2SingleSample -t Plumbing
      else
        echo CloudWorkflow -p SmartSeq2SingleSample ${common_args}
      fi;;
    TargetedSomaticSingleSample)ValidateChip
      continue;;
    ValidateChip)
      echo CloudWorkflow -p ValidateChip ${common_args};;
    VariantCalling)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p VariantCalling -t Plumbing
      else
         echo CloudWorkflow -p VariantCalling ${common_args}
      fi;;
    WholeGenomeGermlineSingleSample)
      echo CloudWorkflow -p WholeGenomeGermlineSingleSample ${common_args};;
    WholeGenomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo CloudWorkflow -p WholeGenomeReprocessing -t Plumbing
      else
        continue
      fi;;
    UltimaGenomicsWholeGenomeGermline)
      echo CloudWorkflow -p UltimaGenomicsWholeGenomeGermline ${common_args};;
    UltimaGenomicsJointGenotyping)
      echo CloudWorkflow -p UltimaGenomicsJointGenotyping ${common_args};;
  esac
}

function main() {
  local -r gittish=${1}
  local -r test_all=${2}
  local -r test=${3}

  local -a changed_pipeline_paths=()
  local -a args=()

  if ${test_all}; then
    changed_pipeline_paths=(${ALL_PIPELINES[@]})
  else
    changed_pipeline_paths=("$(get_pipelines_to_test ${gittish})")
  fi

  for changed_pipeline_path in ${changed_pipeline_paths[*]}; do
    pipeline=$(basename ${changed_pipeline_path} .wdl)
    arg="$(pipeline_to_args ${pipeline} ${test})"
    if [[ -n ${arg} ]]; then
      args+=("${arg}")
    fi
  done

  for arg in "${args[@]}"; do
    echo "${arg}"
  done
}

main ${1} ${2} ${3}