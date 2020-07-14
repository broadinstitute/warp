#!/usr/bin/env bash

set -e

# NOTE: this script will execute from the repository root when called on jenkins

PIPELINE_FOLDER_NAME=$1
TEST_TYPE=${2:-pr}
TEST_TYPE_UC=$(printf $TEST_TYPE | awk '{ print toupper($0) }')
WD=$(pwd)


INPUTS_JSON="/working/test/${PIPELINE_FOLDER_NAME}/${TEST_TYPE}/test_inputs.json"
WDL_FILE="/working/test/${PIPELINE_FOLDER_NAME}/${TEST_TYPE}/test_${PIPELINE_FOLDER_NAME}_${TEST_TYPE_UC}.wdl"
DEPENDENCIES_JSON="/working/test/${PIPELINE_FOLDER_NAME}/${TEST_TYPE}/dependencies.json"

echo "Setting Cromwell environmental variables"

if [ -z ${BROAD_CROMWELL_KEY+x} ]; then
    echo "Error: BROAD_CROMWELL_KEY is not set!"
    exit 1
fi

echo ${BROAD_CROMWELL_KEY} > caas-prod.json
CROMWELL_KEY_FILE="caas-prod.json"

OPTIONS_FILE="/working/test/options.json"
CROMWELL_URL="https://cromwell.caas-prod.broadinstitute.org"
COLLECTION="pipeline-surge"

echo "Running test"
docker run --rm \
  -v "${WD}":/working \
  -w /working \
  --privileged \
  quay.io/broadinstitute/cromwell-tools:v2.1.0 \
  /working/test/test_cromwell_workflow.sh \
    "${CROMWELL_KEY_FILE}" \
    "${CROMWELL_URL}" \
    "${INPUTS_JSON}" \
    "${WDL_FILE}" \
    "${OPTIONS_FILE}" \
    "${DEPENDENCIES_JSON}" \
    "${COLLECTION}"

