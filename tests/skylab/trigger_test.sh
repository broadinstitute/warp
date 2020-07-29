#!/usr/bin/env bash

set -e

# NOTE: this script will execute from the repository root when called on jenkins

PIPELINE_FOLDER_NAME=$1
TEST_TYPE=${2:-pr}
TEST_TYPE_UC=$(printf $TEST_TYPE | awk '{ print toupper($0) }')
WD=$(pwd)

TEST_WDL_FILE="${WD}/tests/skylab/${PIPELINE_FOLDER_NAME}/${TEST_TYPE}/test_${PIPELINE_FOLDER_NAME}_${TEST_TYPE_UC}.wdl"

# Calling release script to flatten dependencies and imports

TMP_DIR=$(mktemp -d -t release-XXXX)
GIT_HASH=$(git rev-parse --short HEAD)

"${WD}/scripts/build_workflow_release.sh" \
    "${TEST_WDL_FILE}" \
    "${TMP_DIR}" \
    "${GIT_HASH}" \
    "${TEST_TYPE}"

DEPENDENCIES_ZIP="dependencies/workflow_dependencies.zip"
FLATTENED_TEST_WDL="test_${PIPELINE_FOLDER_NAME}_${TEST_TYPE_UC}/test_${PIPELINE_FOLDER_NAME}_${TEST_TYPE_UC}.wdl"

# Moving all files into temp dir
cp "${WD}/tests/skylab/test.options.json" "${TMP_DIR}/test.options.json"
OPTIONS_FILE="test.options.json"

cp "${WD}/tests/skylab/${PIPELINE_FOLDER_NAME}/${TEST_TYPE}/test_inputs.json" "${TMP_DIR}/test_inputs.json"
INPUTS_JSON="test_inputs.json"

cp "${WD}/tests/skylab/test_cromwell_workflow.sh" "${TMP_DIR}/test_cromwell_workflow.sh"

echo "Setting Cromwell environmental variables"

if [ -z ${BROAD_CROMWELL_KEY+x} ]; then
    echo "Error: BROAD_CROMWELL_KEY is not set!"
    exit 1
fi

echo ${BROAD_CROMWELL_KEY} > ${TMP_DIR}/caas-prod.json
CROMWELL_KEY_FILE="caas-prod.json"

CROMWELL_URL="https://cromwell.caas-prod.broadinstitute.org"
COLLECTION="pipeline-surge"

echo "Running test"
docker run --rm \
  -v "${TMP_DIR}":/working \
  -w /working \
  --privileged \
  quay.io/broadinstitute/cromwell-tools:v2.1.0 \
  /working/test_cromwell_workflow.sh \
    "${CROMWELL_KEY_FILE}" \
    "${CROMWELL_URL}" \
    "${INPUTS_JSON}" \
    "${FLATTENED_TEST_WDL}" \
    "${OPTIONS_FILE}" \
    "${DEPENDENCIES_ZIP}" \
    "${COLLECTION}"

rm -rf ${TMP_DIR}