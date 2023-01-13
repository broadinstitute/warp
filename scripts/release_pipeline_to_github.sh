#!/usr/bin/env bash

# release_pipeline_to_github.sh takes a WDL and ENVIRONMENT as input
# and pushes a release of the WDL to gitub.
# Along with the release, it also uploads the WDL and dependencies zip
# as assets so the pipeline can be run.
# This script requires that a valid GITHUB_TOKEN oauth token be set
# as an environment variable
set -e

declare -r REPO="broadinstitute/warp"
declare SCRIPT_DIR
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
source ${SCRIPT_DIR}/common.sh

function release_to_github() {
  local -r pipeline=${1} env=${2}

  if [[ ${env} == "${ENV_DEV}" ]]; then
    dev_release_to_github ${pipeline}
  elif [[ ${env} == "${ENV_PROD}" ]]; then
    prod_release_to_github ${pipeline}
  else
    stderr "Only release to dev and prod are enabled."
    exit 1
   fi
}

function dev_release_to_github() {
  local -r pipeline=${1}
  local -r pipelineName=$(basename ${pipeline} .wdl)
  local -r version="develop"
  local -r releaseName=${pipelineName}_${version}
  local -r targetCommitish="develop"
  local -r prerelease=true

  # if there is already a dev release, the tag and release must be deleted
  local -r previousReleaseID=$(curl \
       --silent \
       -X GET \
       -H "Authorization: token ${GITHUB_TOKEN}" \
       "https://api.github.com/repos/${REPO}/releases/tags/${releaseName}" | jq -r .id )

  if [[ ${previousReleaseID} != null ]]; then
      stderr "Deleting tag for ${releaseName}"
      local -r deleteTagResponse=$(curl \
        --silent \
        -X DELETE \
        -H "Authorization: token ${GITHUB_TOKEN}" \
        "https://api.github.com/repos/${REPO}/git/refs/tags/${releaseName}")
      if [[ -n "${deleteTagResponse}" ]]; then
        stderr "Failed to delete the ${releaseName} tag on Github. Printing the response below:"
        stderr $(cat "${deleteTagResponse}")
        exit 1
      fi

      stderr "Deleting ${releaseName} release"
      local -r deleteReleaseResponse=$(curl \
        --silent \
        -X DELETE \
        -H "Authorization: token ${GITHUB_TOKEN}" \
        "https://api.github.com/repos/${REPO}/releases/${previousReleaseID}")
      if [[ -n "${deleteReleaseResponse}" ]]; then
        stderr "Failed to delete the ${releaseName} release on Github. The tag was deleted, but the release and artifacts remain.
         To proceed, manually delete the release with the id ${previousReleaseID} on github. Then re-run the release script from the
         dsde-pipelines root directory by running './scripts/release_pipeline_to_github.sh -p ${pipeline} -e dev'.
         Printing the response below:"
        stderr $(cat "${deleteReleaseResponse}")
        exit 1
      fi
  fi

  build_and_release_to_github ${pipeline} ${version} ${releaseName} ${targetCommitish} ${prerelease}
}

function prod_release_to_github() {
  local -r pipeline=${1}
  local -r pipelineName=$(basename ${pipeline} .wdl)
  local -r changelog=$(dirname ${pipeline})/${pipelineName}.changelog.md
  local -r version=$(get_version_from_changelog ${changelog})
  local -r releaseName=${pipelineName}_v${version}
  local -r targetCommitish="master"
  local -r prerelease=false

  local -r previousReleaseID=$(curl \
      --silent \
      -X GET \
      -H "Authorization: token ${GITHUB_TOKEN}" \
      "https://api.github.com/repos/${REPO}/releases/tags/${releaseName}" | jq -r .id )
  if [[ ${previousReleaseID} != null ]]; then
      stderr "There is a previous release of ${releaseName} on Github with the release id ${previousReleaseID}."
      exit 0
  fi

  build_and_release_to_github ${pipeline} ${version} ${releaseName} ${targetCommitish} ${prerelease}
}

function build_and_release_to_github() {
  local -r pipeline=${1}
  local -r version=${2}
  local -r releaseName=${3}
  local -r targetCommitish=${4}
  local -r prerelease=${5}
  local -r pipelineName=$(basename ${pipeline} .wdl)
  local -r changelog=$(dirname ${pipeline})/${pipelineName}.changelog.md

  stderr "Building artifacts for ${releaseName}"
  ${SCRIPT_DIR}/build_pipeline_release.sh -w ${pipeline} -e prod -v ${version} -o ${localReleaseDir}

  stderr "Building release notes for ${releaseName}"
  local previousEntryStart
  previousEntryStart=$(grep -n '^# [0-9]*\.[0-9]*' ${changelog} | sed -n 2p | cut -d ':' -f 1)
  if [[ -z ${previousEntryStart} ]]; then
    previousEntryStart=$(wc -l ${changelog} | awk '{print $1}')
  else
    previousEntryStart=$((previousEntryStart - 1))
  fi
  local -r body="$(head -n ${previousEntryStart} "${changelog}" | sed 's/$/\\r\\n/g' | tr -d '\n')"

  local -r payload="{\"tag_name\": \"${releaseName}\",\"name\": \"${releaseName}\",\"target_commitish\": \"${targetCommitish}\",\"body\": \"${body}\",\"prerelease\": ${prerelease}}"

  stderr "Making the ${releaseName} release on Github"

  local -r releaseId=$(upload_to_github_as_draft \
    ${pipelineName} \
    ${version} \
    ${releaseName} \
    "${payload}")

  if [[ -z ${releaseId} ]]; then
    stderr "Failed to make a draft release to Github"
    exit 1
  else
    stderr "Draft release ${releaseName} was successfully uploaded to Github. Now publishing..."
    publish_to_github ${releaseId}
    stderr "Done making the ${releaseName} release on Github"
    rm -rf ${localReleaseDir}
  fi
}

function publish_to_github() {
  local -r releaseId=${1}
  local publishResponse=${localReleaseDir}/publishResponse

  local -r publishResponseCode=$(curl \
    --silent \
    --output ${publishResponse} \
    --write-out "%{http_code}" \
    -X PATCH \
    -H "Authorization: token ${GITHUB_TOKEN}" \
    -d '{"draft":false}' \
    "https://api.github.com/repos/${REPO}/releases/${releaseId}")

  if [[ ${publishResponseCode} -ne 200 ]]; then
    stderr "Failed to publish the release on GitHub (${publishResponseCode}). Printing the response below:"
    stderr $(cat "${publishResponse}")
    exit 1
  fi
}

function upload_to_github_as_draft() {
  local -r pipelineName=${1}
  local -r version=${2}
  local -r releaseName=${3}
  local -r payload="${4}"

  local -r releaseResponse=${localReleaseDir}/releaseResponse

  local -r releaseResponseCode=$(curl \
    --silent \
    --output ${releaseResponse} \
    --write-out "%{http_code}" \
    -X POST \
    -H "Authorization: token ${GITHUB_TOKEN}" \
    -d "$(echo "${payload}" | jq '. + {"draft":true}')" \
    "https://api.github.com/repos/${REPO}/releases")

  if [[ ${releaseResponseCode} -ne 201 ]]; then
    stderr "Failed to create the draft release on Github. Printing the response below:"
    stderr "${releaseResponseCode}"
    stderr $(cat "${releaseResponse}")
    exit 1
  fi

  local -r releaseId=$(jq -r .id ${releaseResponse})

  upload_artifact_to_github \
    ${releaseId} \
    ${pipelineName} \
    ${version} \
    ${releaseName} \
    "wdl" \
    "text/plain"

  upload_artifact_to_github \
    ${releaseId} \
    ${pipelineName} \
    ${version} \
    ${releaseName} \
    "options.json" \
    "application/json"

  local -r dependenciesZip=${localReleaseDir}/${pipelineName}/${pipelineName}_${version}.zip

  if [[ -f ${dependenciesZip} ]]; then
    upload_artifact_to_github \
      ${releaseId} \
      ${pipelineName} \
      ${version} \
      ${releaseName} \
      "zip" \
      "application/zip"
  fi

  echo ${releaseId}
}

function upload_artifact_to_github() {
  local -r releaseId=${1}
  local -r pipelineName=${2}
  local -r version=${3}
  local -r releaseName=${4}
  local -r artifact=${5}  # 'zip', 'options.json', or 'wdl'
  local -r contentType=${6} # 'text/plain' or 'application/zip'


  stderr "Uploading ${releaseName}.${artifact} to Github"

  local -r response=${localReleaseDir}/${artifact}Response
  local -r responseCode=$(curl \
      --silent \
      --output ${response} \
      --write-out "%{http_code}" \
      -X POST \
      -H "Authorization: token ${GITHUB_TOKEN}" \
      -H "Content-Type: ${contentType}" \
      --data-binary @"${localReleaseDir}/${pipelineName}/${pipelineName}_${version}.${artifact}" \
      "https://uploads.github.com/repos/${REPO}/releases/${releaseId}/assets?name=${releaseName}.${artifact}")

  if [[ ${responseCode} -ne 201 ]]; then
    stderr "Failed to upload ${artifact} to Github. Printing the response below:"
    stderr $(cat "${response}")
    exit 1
  fi
}

function cleanup_failed_release() {
  if [[ -f ${localReleaseDir}/releaseResponse ]]; then
    local -r releaseId=$(jq -r .id ${localReleaseDir}/releaseResponse)
    if [[ ${releaseId} != "null" ]]; then
      stderr "Cleaning up failed release on Github"
      curl \
        --silent \
        -X DELETE \
        -H "Authorization: token ${GITHUB_TOKEN}" \
        "https://api.github.com/repos/${REPO}/releases/${releaseId}" >&2

      local -r tagName=$(jq -r .tag_name ${localReleaseDir}/releaseResponse)
      stderr "Cleaning up the tag for failed release on Github"
      curl \
        --silent \
        -X DELETE \
        -H "Authorization: token ${GITHUB_TOKEN}" \
        "https://api.github.com/repos/${REPO}/git/refs/tags/${tagName}" >&2
    else
      stderr "Releasing failed before anything could be uploaded."
    fi
    rm -rf ${localReleaseDir}
    exit 1
  fi
}

declare PIPELINE_TO_RELEASE=""
declare ENV=${ENV_DEV}

while getopts "p:e:" opt; do
  case ${opt} in
    p)
      if [[ ! -f ${OPTARG} ]]; then
        stderr "Error: ${OPTARG} does not exist!"
        exit 1
      fi
      PIPELINE_TO_RELEASE=${OPTARG}
      ;;
    e)
      if [[ ${ENVIRONMENTS[*]} =~ ${OPTARG} ]]; then
        ENV=${OPTARG}
      else
        stderr "usage: -e [${ENVIRONMENTS[*]}]"
        exit 1
      fi
      ;;
    *)
      >&2 stderr "invalid option"
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))

if [[ -z ${GITHUB_TOKEN} ]]; then
  stderr "GITHUB_TOKEN must be set as an environment variable"
  exit 1
fi

declare localReleaseDir
localReleaseDir=$(mktemp -d)

trap cleanup_failed_release EXIT

release_to_github ${PIPELINE_TO_RELEASE} ${ENV}
