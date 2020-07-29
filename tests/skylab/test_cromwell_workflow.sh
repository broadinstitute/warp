#!/usr/bin/env bash

set -e

CROMWELL_KEY_FILE=$1
CROMWELL_URL=$2
INPUTS_JSON=$3
WDL_FILE=$4
OPTIONS_FILE=$5
DEPENDENCIES_ZIP=$6
COLLECTION=$7

echo "Starting workflow."

WORKFLOW_HASH=$(cromwell-tools submit \
  --service-account-key "${CROMWELL_KEY_FILE}" \
  --url "${CROMWELL_URL}" \
  --inputs-files "${INPUTS_JSON}" \
  --wdl-file "${WDL_FILE}" \
  --options-file "${OPTIONS_FILE}" \
  --collection-name "${COLLECTION}" \
  --deps-file "${DEPENDENCIES_ZIP}" \
)

# Get workflow id from cromwell-tools response: {"id": "XXXXX", "status": "Submitted"}
WORKFLOW_ID=$(echo ${WORKFLOW_HASH} | python3 -c "import json,sys;obj=json.load(sys.stdin);print(obj['id']);")
echo ${WORKFLOW_ID}

# Wait before polling because there is a delay in Cromwell when updating the workflow status
# Note: This can be removed if cromwell-tools retries this API call internally:
# https://github.com/broadinstitute/cromwell-tools/issues/62
sleep 10

cromwell-tools --version

echo "Waiting for workflow ${WORKFLOW_ID} to complete..."
cromwell-tools wait "${WORKFLOW_ID}" \
  --service-account-key "${CROMWELL_KEY_FILE}" \
  --url "${CROMWELL_URL}" \
  --timeout-minutes 180 \
  --poll-interval-seconds 30

echo "Checking workflow completion status."
STATUS_RESULT=$(cromwell-tools status \
  --service-account-key "${CROMWELL_KEY_FILE}" \
  --url "${CROMWELL_URL}" \
  --uuid "${WORKFLOW_ID}" \
)

# get the status alone; cromwell tools is a bit more verbose
STATUS=$(echo ${STATUS_RESULT} | python3 -c "import json,sys;obj=json.load(sys.stdin);print(obj['status']);")

if [ "${STATUS}" == "Succeeded" ]; then
  echo "Test Successful. Workflow outputs matched expected values."
  exit 0
else
  echo "Test Failed. Workflow ${WORKFLOW_HASH} did not match expected outputs. See logs for listed workflow for more details."
  exit 1
fi
