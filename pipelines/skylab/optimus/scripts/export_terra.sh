#!/bin/bash

declare -r namespace=${WORKSPACE_NAMESPACE}
declare -r name=${WORKSPACE_NAME}

python3 workspace_controller.py get_participant_table \
--workspace-name $name \
--participant-table-name participant \
--output-name participant.tsv

python3 workspace_controller.py create_participant_lane \
--input-name participant.tsv \
--output-prefix participant_lane

python3 workspace_controller.py upload_participant \
--workspace-name $name \
--input-prefix participant_lane

echo "Data uploaded to ${namespace}/${name}"
