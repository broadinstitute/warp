#!/usr/bin/env bash

# This script takes a commit/tag/ref as an argument and scans all WDLs in the repo
# to see if they require changes as outlined by our contributor guidelines.
# It is run during PRs and fails if pipelines with changed dependencies have
# not been updated.
# You can optionally supply a boolean using the -i parameter to ignore version
# numbers and only check the changelogs. For instance '-i true' will result
# in only changelogs being checked.

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/common.sh

function check_if_file_changed() {
  local -r file=${1} commit=${2}
  # the result is flipped because git diff-index returns 0 if the file is not changed
  ! git diff-index --quiet ${commit} ${file}
}

function get_version_from_workflow() {
  local -r workflow=${1}
  echo $(grep -m1 '^\s.*String pipeline_version\s.*=\s"[0-9]*\.[0-9]*' ${workflow} | cut -f 2 -d '=' | xargs)
}

function validate_version_changed() {
  local -r workflow=${1} commit=${2}
  local -r linesChanged=$(git diff ${commit} -- \
    ${workflow} |
    grep '^[+]\s.*String pipeline_version\s.*=\s"[0-9]*\.[0-9]*' |
    wc -l)
  test ${linesChanged} -eq 1
}

function validate_wdl_changes() {
  local -r workflow=${1} commitToCompare=${2}
  local passesValidation=0

  # First check if the WDL has been changed
  if $(check_if_file_changed ${workflow} ${commitToCompare}); then
    changelog=$(dirname ${workflow})/$(basename ${workflow} .wdl).changelog.md

    # Check if pipeline version has been updated
    if ! $(validate_version_changed ${workflow} ${commitToCompare}); then
      stderr "$(basename ${workflow}) has not had its version updated"
      passesValidation=1
    fi

    # Check if the pipeline version matches the changelog version
    local -r changelogVersion=$(get_version_from_changelog ${changelog})
    local -r pipelineVersion=$(get_version_from_workflow ${workflow})
    if [[ ${pipelineVersion} != ${changelogVersion} ]]; then
      stderr "$(basename ${workflow}) and $(basename ${changelog}) do not have matching versions"
      passesValidation=1
    fi
  else
    stderr "$(basename ${workflow}) has not been changed and needs updating"
    passesValidation=1
  fi

  return ${passesValidation}
}

function validate_changelog_changes() {
  local -r changelog=${1} commitToCompare=${2}
  local validChangelog=0

  if ! $(check_if_file_changed ${changelog} ${commitToCompare}); then
    stderr "$(basename ${changelog}) has not been changed and needs to be updated"
    validChangelog=1
  fi
  return ${validChangelog}
}


function validate_versions_and_changelogs() {
  local commitToCompare=${1}
  local -r -a modifiedPipelines=($(get_modified_pipelines ${commitToCompare}))
  local allPipelinesValid=true

  stderr "Comparing versions and changelogs for pipelines that differ from the versions on '${commitToCompare}':"
  for modifiedPipeline in ${modifiedPipelines[@]}; do
    if ! $(validate_wdl_changes ${modifiedPipeline} ${commitToCompare}); then
      allPipelinesValid=false
    fi
  done

  if ${allPipelinesValid}; then
    stderr "All WDLs and changelog files appear to be valid for this release."
    return 0
  else
    stderr "Some WDLs or changelog files need updating. See output for details."
    return 1
  fi

}

function validate_changelogs_only() {
  local commitToCompare=${1}
  local -r -a modifiedPipelines=($(get_modified_pipelines ${commitToCompare}))
  local allPipelinesValid=true

  stderr "Comparing changelogs for pipelines that differ from the versions on '${commitToCompare}':"
  for modifiedPipeline in ${modifiedPipelines[@]}; do
    changelog=$(dirname ${modifiedPipeline})/$(basename ${modifiedPipeline} .wdl).changelog.md
    if ! $(validate_changelog_changes ${changelog} ${commitToCompare}); then
      allPipelinesValid=false
    fi
  done

  if ${allPipelinesValid}; then
    stderr "All changelog files are valid for this release."
    return 0
  else
    stderr "Some changelog files need updating. See output for details."
    return 1
  fi
}

function validate_pipelines() {
  local commitToCompare=${1}

  if ${IGNORE_VERSION}; then
    return $(validate_changelogs_only ${commitToCompare})
  else
    return $(validate_versions_and_changelogs ${commitToCompare})
  fi
}

declare GIT_COMMITISH=""
declare IGNORE_VERSION=false

while getopts "g:i:" opt; do
  case ${opt} in
    g)
      GIT_COMMITISH=${OPTARG}
      ;;
    i)
      IGNORE_VERSION=${OPTARG}
      ;;
    *)
      >&2 echo "invalid option"
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))

if [[ -z ${GIT_COMMITISH} ]]; then
  stderr "You must supply a git commitish using the -g parameter."
fi

exit $(validate_pipelines ${GIT_COMMITISH})
