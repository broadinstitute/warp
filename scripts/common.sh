#!/usr/bin/env bash

declare -r REPO_ROOT=$(git rev-parse --show-toplevel)
declare -r -a ENVIRONMENTS=(dev rc prod)
declare -r ENV_DEV=${ENVIRONMENTS[0]}
declare -r ENV_RC=${ENVIRONMENTS[1]}
declare -r ENV_PROD=${ENVIRONMENTS[2]}

function stderr() {
  >&2 echo "$@"
}

function get_version_from_changelog() {
  local -r changelog=${1}
  echo $(grep -m1 '^#' ${changelog} | cut -f 2 -d ' ')
}

function get_dependencies_for_wdl() {
  local -r wdl=${1}
  local -r wdlDir=$(dirname "${wdl}")
  local -a wdlImports=()

  while IFS= read -r importPath; do
    # Extract import path and remove quotes
    cleanPath=$(echo "${importPath}" | awk '{print $2}' | tr -d '"')

    # Skip http imports
    if [[ "${cleanPath}" =~ ^http ]]; then
      continue
    fi

    # Resolve to absolute path relative to the WDL’s directory
    local resolvedPath
    resolvedPath=$(realpath -m "${wdlDir}/${cleanPath}")

    # Only add if the file exists
    if [[ -f "${resolvedPath}" ]]; then
      wdlImports+=("${resolvedPath}")
    else
      stderr "Warning: could not resolve import ${cleanPath} from ${wdl}"
    fi
  done < <(grep '^import ".*$' "${wdl}")

  local -a subWorkflowImports=()
  for import in "${wdlImports[@]}"; do
    subWorkflowImports+=($(get_dependencies_for_wdl "${import}"))
  done

  echo "${wdlImports[@]}" "${subWorkflowImports[@]}"
}


function get_pipeline_dependencies() {
  local -r pipeline=${1}
  local -a basenameDependencies=()
  local -r dependencies=($(get_dependencies_for_wdl ${pipeline}))
  for dependency in ${dependencies[@]}; do
    basenameDependencies+=($(basename ${dependency}))
  done

  echo ${basenameDependencies[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '
}

function get_versioned_pipelines() {
  local -r -a pipelines=($(find ${REPO_ROOT}/pipelines -type f -name '*.wdl'))
  echo ${pipelines[@]}
}

function get_modified_pipelines() {
  local -r commitToCompare=${1}
  local -r -a pipelines=($(get_versioned_pipelines))
  local -r -a changedWdls=($(git diff --name-only HEAD ${commitToCompare} | grep -E '.*\.wdl' | xargs -n1 basename))
  local -a modifiedPipelines=()
  for pipeline in ${pipelines[@]}; do
    if [[ " ${changedWdls[@]}" =~ " $(basename ${pipeline})" ]]; then
      modifiedPipelines+=(${pipeline})
    else
      dependencies=($(get_pipeline_dependencies ${pipeline}))
      for changedWdl in ${changedWdls[@]}; do
        if [[ " ${dependencies[@]}" =~ " ${changedWdl}" ]]; then
          modifiedPipelines+=(${pipeline})
          break
        fi
      done
    fi
  done
  echo ${modifiedPipelines[@]}
}

function get_pipelines_to_test() {
  local -r commitToCompare=${1}
  local -r -a pipelines=($(get_versioned_pipelines))
  local -r -a changedFiles=($(git diff --name-only HEAD ${commitToCompare}))
  local -a pipelinesToTest=()

  if [[ "${changedFiles[@]}" == *"verification/"*  || \
        "${changedFiles[@]}" == *"tests/broad"* ]]; then
    pipelinesToTest+=(${pipelines[@]})
  else
    local -a modifiedPipelines=($(get_modified_pipelines ${commitToCompare}))
    for pipeline in ${pipelines[@]}; do
      # Name comparison includes a leading space to avoid discovering pipelines named with substrings of other pipeline names, i.e. "Arrays" and  "MutilSampleArrays"
      if [[ " ${modifiedPipelines[@]}" =~ " ${pipeline}" ]]; then
        pipelinesToTest+=(${pipeline})
      else
        pipelineDirectory=$(dirname ${pipeline})
        pipelineHomeDir=${pipelineDirectory/${REPO_ROOT}\//}
        for changedFile in ${changedFiles[@]}; do
          if [[ "$(dirname ${changedFile})" == "${pipelineHomeDir}"* ]]; then
            pipelinesToTest+=(${pipeline})
            break
          fi
        done
      fi
    done
  fi

  echo ${pipelinesToTest[@]}
}
