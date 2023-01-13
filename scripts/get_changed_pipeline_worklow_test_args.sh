#!/usr/bin/env bash

declare SCRIPT_DIR
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
source ${SCRIPT_DIR}/common.sh

declare ALL_PIPELINES=($(get_versioned_pipelines))

function pipeline_to_args() {
  local -r pipeline=${1}
  local -r test=${2}

  local -r common_args="${test}"

  case ${pipeline} in
    AnnotationFiltration)
      continue;;
    Arrays)
      echo Arrays ${common_args};;
    BroadInternalRNAWithUMIs)
      echo BroadInternalRNAWithUMIs ${common_args};;
    BroadInternalUltimaGenomics)
      echo BroadInternalUltimaGenomics ${common_args};;
    # CEMBA)
    #   echo CEMBA ${common_args};;
    CheckFingerprint)
      echo CheckFingerprint ${common_args};;
    CramToUnmappedBams)
      echo CramToUnmappedBams ${common_args};;
    ExternalExomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo ExternalExomeReprocessing Plumbing
      else
        continue
      fi;;
    ExomeGermlineSingleSample)
      echo ExomeGermlineSingleSample ${common_args};;
    ExomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo ExomeReprocessing Plumbing
      else
        continue
      fi;;
    ExternalWholeGenomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo ExternalWholeGenomeReprocessing Plumbing
      else
        continue
      fi;;
    GDCWholeGenomeSomaticSingleSample)
      echo GDCWholeGenomeSomaticSingleSample ${common_args};;
    IlluminaGenotypingArray)
      echo IlluminaGenotypingArray ${common_args};;
    Imputation)
      echo Imputation ${common_args};;
    JointGenotyping)
      echo JointGenotyping ${common_args};;
    JointGenotypingByChromosomePartOne)
      continue;;
    JointGenotypingByChromosomePartTwo)
      continue;;
    MultiSampleArrays)
      echo MultiSampleArrays ${common_args};;
    MultiSampleSmartSeq2)
      if [[ "${test}" == "Scientific" ]]; then
        echo MultiSampleSmartSeq2 Plumbing
      else
        echo MultiSampleSmartSeq2 ${common_args}
      fi;;
    MultiSampleSmartSeq2SingleNucleus)
      if [[ "${test}" == "Scientific" ]]; then
        echo MultiSampleSmartSeq2SingleNucleus Plumbing
      else
        echo MultiSampleSmartSeq2SingleNucleus ${common_args}
      fi;;
    Optimus)
      echo Optimus ${common_args};;
    ReblockGVCF)
      echo ReblockGvcf  ${common_args};;
    RNAWithUMIsPipeline)
      echo RNAWithUMIsPipeline ${common_args};;
    scATAC)
      if [[ "${test}" == "Scientific" ]]; then
        echo scATAC Plumbing
      else
        echo scATAC ${common_args}
      fi;;
    SmartSeq2SingleSample)
      if [[ "${test}" == "Scientific" ]]; then
        echo SmartSeq2SingleSample Plumbing
      else
        echo SmartSeq2SingleSample ${common_args}
      fi;;
    TargetedSomaticSingleSample)ValidateChip
      continue;;
    ValidateChip)
      echo ValidateChip ${common_args};;
    VariantCalling)
      if [[ "${test}" == "Scientific" ]]; then
        echo VariantCalling Plumbing
      else
         echo VariantCalling ${common_args}
      fi;;
    WholeGenomeGermlineSingleSample)
      echo WholeGenomeGermlineSingleSample ${common_args};;
    WholeGenomeReprocessing)
      if [[ "${test}" == "Scientific" ]]; then
        echo WholeGenomeReprocessing Plumbing
      else
        continue
      fi;;
    UltimaGenomicsWholeGenomeCramOnly)
      echo UltimaGenomicsWholeGenomeCramOnly ${common_args};;
    UltimaGenomicsWholeGenomeGermline)
      echo UltimaGenomicsWholeGenomeGermline ${common_args};;
    UltimaGenomicsJointGenotyping)
      echo UltimaGenomicsJointGenotyping ${common_args};;
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