#!/usr/bin/env bash

set -e

OLD_DIR=${PWD}
SCRIPT_DIR=$(dirname $(readlink -f "$0"))

cd ${SCRIPT_DIR}

function finish {
  cd ${OLD_DIR}
}

trap finish EXIT

DSDE_PIPELINES_DIR=$(git rev-parse --show-toplevel)

## PROJECT-SPECIFIC CONSTANTS

# Don't put a '/' at the end of the OUTPUT_CLOUD_DIR
declare -r OUTPUT_CLOUD_DIR=""
declare -r CALLSET_NAME=""
declare -r SAMPLE_NAME_MAP="/change_me"
declare -r DATA_TYPE="{{DATA_TYPE}}" # DATA_TYPE should be 'wgs' or 'exome'

## END PROJECT-SPECIFIC CONSTANTS. DO NOT EDIT CONSTANTS BELOW THIS LINE

declare -r CROMWELL_URL="https://cromwell-jg.gotc-prod.broadinstitute.org/api/workflows/v1"
declare -r GCS_EXECUTION_BUCKET="gs://broad-gotc-prod-execution"
declare -r PART_ONE_FOFN_DIR="$(dirname ${SAMPLE_NAME_MAP})/part_one_fofn"
declare -r COMPUTE_PROJECT_PREFIX="{{COMPUTE_PROJECT_PREFIX}}"
declare -r PROJECT_NUMBER_RANGE=({{PROJECT_NUMBER_RANGE_START}} {{PROJECT_NUMBER_RANGE_END}})
declare -r PART_ONE_OUTPUTS_SUBDIR="part_one_outputs"


# Creates a working copy of this script and uses sed to replace variables in the working copy
# with the appropriate values depending on the data type.
# For instance, if this script is called with the argument 'WGS', all instances of '{{DATA_TYPE}}' will be changed to "wgs",
# all instances of '{{COMPUTE_PROJECT_PREFIX}}' will be replaced with "broad-wgs-jg-prod", and so on
function do_setup() {
  case ${1} in
  WGS)
    dataType="wgs"
    computeProjectPrefix="broad-wgs-jg-prod"
    projectNumberRangeStart="1"
    projectNumberRangeEnd="2"
    ;;
  Exome)
    dataType="exome"
    computeProjectPrefix="broad-exomes-prod"
    projectNumberRangeStart="1"
    projectNumberRangeEnd="10"
    ;;
    *)
    >&2 echo "${1} is not a valid data type. Must be 'WGS' or 'Exome'"
    exit 1
    ;;
esac

  >&2 echo "Setting up shop in ${OLD_DIR}"
  ln -fs ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/${dataType}/JointGenotypingByChromosomePartOne.template.input.json ${OLD_DIR}/JointGenotypingByChromosomePartOne.input.json
  ln -fs ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/${dataType}/JointGenotypingByChromosomePartTwo.template.input.json ${OLD_DIR}/JointGenotypingByChromosomePartTwo.input.json

  cp \
    ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/by_chromosome_client.sh \
    ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/by_chromosome_client.working.sh

  sed -i.bak 's|{{DATA_TYPE}}|'${dataType}'|g' ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/by_chromosome_client.working.sh
  sed -i.bak 's|{{COMPUTE_PROJECT_PREFIX}}|'${computeProjectPrefix}'|g' ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/by_chromosome_client.working.sh
  sed -i.bak 's|{{PROJECT_NUMBER_RANGE_START}}|'${projectNumberRangeStart}'|g; s|{{PROJECT_NUMBER_RANGE_END}}|'${projectNumberRangeEnd}'|g' ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/by_chromosome_client.working.sh
  ln -fs ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/by_chromosome_client.working.sh ${OLD_DIR}/by_chromosome_client.sh

  rm ${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/*.bak
}

if [[ ${1} == "setup" ]]; then
  do_setup ${2}
  exit 0
fi

JOINT_GENOTYPING_PART_ONE_WDL=${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/JointGenotypingByChromosomePartOne.wdl
JOINT_GENOTYPING_PART_TWO_WDL=${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/JointGenotypingByChromosomePartTwo.wdl

JOINT_GENOTYPING_PART_ONE_INPUTS_TEMPLATE=${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/${DATA_TYPE}/JointGenotypingByChromosomePartOne.template.input.json
JOINT_GENOTYPING_PART_TWO_INPUTS_TEMPLATE=${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/${DATA_TYPE}/JointGenotypingByChromosomePartTwo.template.input.json

CROMWELL_OPTIONS_TEMPLATE=${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/JointGenotypingByChromosome.template.options.json

WORKFLOWS_CSV=${SCRIPT_DIR}/workflows.csv
CHROMOSOME_SCATTER_CSV=${DSDE_PIPELINES_DIR}/pipelines/broad/dna_seq/germline/joint_genotyping/by_chromosome/${DATA_TYPE}/chromosome_scatter_counts.csv
RELEASE_DIR=''

KEYS_TO_INCLUDE="includeKey=executionStatus&includeKey=backendStatus&includeKey=status&includeKey=workflowName&includeKey=subWorkflowMetadata&includeKey=subWorkflowId&includeKey=id"
COST_KEYS_TO_EXCLUDE="excludeKey=submittedFiles&excludeKey=inputs&excludeKey=outputs&excludeKey=commandLine&excludeKey=backendLabels&excludeKey=callRoot&excludeKey=stderr&excludeKey=stdout&excludeKey=dockerImageUsed&excludeKey=labels&excludeKey=jobId&excludeKey=backendLogs&excludeKey=hashes&excludeKey=input&excludeKey=endpointUrl&excludeKey=monitoringScript&excludeKey=executionBucket&excludeKey=instanceName&excludeKey=zone&excludeKey=cpuMin&excludeKey=failOnStderr&excludeKey=memoryMin&excludeKey=continueOnReturnCode&excludeKey=maxRetries&excludeKey=noAddress&excludeKey=allowResultReuse"
EXPAND_SUB_WORKFLOWS_ARG="expandSubWorkflows=true"
KEYS_TO_EXCLUDE="excludeKey=submittedFiles"

# Do a POST request
function do_post() {
  curl --silent -X POST "$@"
}

# Do a GET request
function do_get() {
  curl --silent -X GET "$@"
}

# Get a workflow ID from the log
# `assert_id_exists` should be called after this
# They can't be combined because exiting subshells is hard
function get_workflow_id_from_log() {
  local -r chromosome=${1}
  if ! grep "${chromosome}," ${WORKFLOWS_CSV} > /dev/null; then
    echo ""
  fi
  local -r id=$(grep "^${chromosome}," ${WORKFLOWS_CSV} | cut -d ',' -f 2)
  echo ${id}
}

# Assert that the ID just gotten from `get_workflow_id_from_log` is real
# Exit if its not real
function assert_id_exists() {
  local -r chromosome=${1} id=${2}
  if [[ -z "${id}" ]]; then
    >&2 echo "There is no workflow ID for ${chromosome}"
    exit 1
  fi
}

# Release a workflow to a temporary directory
# This bundles the workflow dependencies necessary for submission
function do_release() {
  local -r wdl=${1}
  RELEASE_DIR=$(mktemp -d)

  ${DSDE_PIPELINES_DIR}/scripts/build_workflow_release.sh \
    ${wdl} \
    ${RELEASE_DIR} \
    $(git rev-parse --short HEAD) \
    pharma5 > /dev/null

  echo Released to temporary directory ${RELEASE_DIR}
}

# Submit a chromosome for processing
function do_submit_part_one() {
  local -r chromosome=${1} projectNum=${2}

  local -r tmpInputs=$(mktemp)
  local -r tmpOptions=$(mktemp)

  do_release ${JOINT_GENOTYPING_PART_ONE_WDL}
  local -r wdl=${RELEASE_DIR}/JointGenotypingByChromosomePartOne/JointGenotypingByChromosomePartOne.wdl

  local -r chromosomeScatter=$(grep "${chromosome}," ${CHROMOSOME_SCATTER_CSV} | cut -d ',' -f 2)
  sed 's|{CHROMOSOME}|'${chromosome}'|g; s|{CHROMOSOME_SCATTER}|'${chromosomeScatter}'|g; s|{CALLSET_NAME}|'${CALLSET_NAME}'|g; s|{SAMPLE_NAME_MAP}|'${SAMPLE_NAME_MAP}'|g' ${JOINT_GENOTYPING_PART_ONE_INPUTS_TEMPLATE} > ${tmpInputs}
  sed 's|{RANDOM_PROJECT}|'${projectNum}'|g; s|{COMPUTE_PROJECT_PREFIX}|'${COMPUTE_PROJECT_PREFIX}'|g; s|{GCS_EXECUTION_BUCKET}|'${GCS_EXECUTION_BUCKET}'|g' ${CROMWELL_OPTIONS_TEMPLATE} > ${tmpOptions}

  echo "Submitting ${chromosome} to cromwell"
  submit_job ${wdl} ${chromosome} ${tmpInputs} ${tmpOptions} ${projectNum}

  rm ${tmpInputs}
  rm ${tmpOptions}
}

# Submit the results of all the part one chromosomes for processing
function do_submit_part_two() {
  local -r projectNum=${1}

  # Collect all the vcfs from part one and provide them as input to part two
  write_part_two_inputs_to_cloud

  local -r tmpInputs=$(mktemp)
  local -r tmpOptions=$(mktemp)

  do_release ${JOINT_GENOTYPING_PART_TWO_WDL}
  local -r wdl=${RELEASE_DIR}/JointGenotypingByChromosomePartTwo/JointGenotypingByChromosomePartTwo.wdl

  sed 's|{PART_ONE_FOFN_DIR}|'${PART_ONE_FOFN_DIR}'|g; s|{SAMPLE_NAME_MAP}|'${SAMPLE_NAME_MAP}'|g; s|{CALLSET_NAME}|'${CALLSET_NAME}'|g' ${JOINT_GENOTYPING_PART_TWO_INPUTS_TEMPLATE} > ${tmpInputs}

  sed 's|{RANDOM_PROJECT}|'${projectNum}'|g; s|{COMPUTE_PROJECT_PREFIX}|'${COMPUTE_PROJECT_PREFIX}'|g; s|{GCS_EXECUTION_BUCKET}|'${GCS_EXECUTION_BUCKET}'|g' ${CROMWELL_OPTIONS_TEMPLATE} > ${tmpOptions}

  echo "Submitting JointGenotypingByChromosomePartTwo to cromwell"
  submit_job ${wdl} "part-two" ${tmpInputs} ${tmpOptions} ${projectNum}

  rm ${tmpInputs} ${tmpOptions}
}

# Submit a WDL to Cromwell.
function submit_job() {
  local -r wdl=${1} label=${2} inputs=${3} options=${4} projectNum=${5}

  if [[ ! -z $(grep "${label}," ${WORKFLOWS_CSV} | cut -d ',' -f 2) ]]; then
    >&2 echo "An there is already an entry for ${label}."
    >&2 echo "Please abort the current workflow or mark it as completed."
    exit 1
  fi

  local -r response=$(do_post \
    -F workflowSource=@${wdl} \
    -F workflowInputs=@${inputs} \
    -F workflowOptions=@${options} \
    -F workflowDependencies=@${RELEASE_DIR}/dependencies/workflow_dependencies.zip \
    -F labels='{"chromosome":"'${label}'"}' \
    ${CROMWELL_URL} 2>/dev/null)

  local -r status=$(echo ${response} | jq -r .status)
  local -r id=$(echo ${response} | jq -r .id)
  if grep "${label}," ${WORKFLOWS_CSV} >/dev/null; then
    echo "Updating entry for ${label} entry"
    sed -i.bak -E 's/^'${label}',([a-f0-9-]*),[0-9]+,([|a-f0-9-]*)/'${label}','${id}','${projectNum}',\2\1/g' ${WORKFLOWS_CSV}
    rm ${WORKFLOWS_CSV}.bak
  else
    echo "Adding ${label} entry to the log"
    echo -e "${label},${id},${projectNum}," >> ${WORKFLOWS_CSV}

  fi
  echo "Workflow ${id} is now ${status} in ${COMPUTE_PROJECT_PREFIX}${projectNum}"
}

# Write the outputs of part one to a cloud dir for part two's consumption
function write_part_two_inputs_to_cloud() {
  local -r sitesOnlyVcfs=$(mktemp)
  local -r sitesOnlyVcfIndices=$(mktemp)
  local -r annotationVcfs=$(mktemp)
  local -r annotationVcfIndices=$(mktemp)
  local -r hardFilteredVcfs=$(mktemp)
  local -r hardFilteredVcfIndices=$(mktemp)
  local -r fingerprintingVcfs=$(mktemp)
  local -r fingerprintingVcfIndices=$(mktemp)
  while IFS=, read chromosome workflowId project failedIds; do
    if [[ ${chromosome} == chr* ]]; then
      >&2 echo "Collecting outputs for ${chromosome}"
      outputs=$(mktemp)
      gsutil ls "${OUTPUT_CLOUD_DIR}/${PART_ONE_OUTPUTS_SUBDIR}/${chromosome}/" > ${outputs}
      grep 'sites_only.vcf.gz$' ${outputs} >> ${sitesOnlyVcfs}
      grep 'sites_only.vcf.gz.tbi$' ${outputs} >> ${sitesOnlyVcfIndices}
      grep 'annotationDB.vcf.gz$' ${outputs} >> ${annotationVcfs}
      grep 'annotationDB.vcf.gz.tbi$' ${outputs} >> ${annotationVcfIndices}
      grep 'hard_filtered_with_genotypes.vcf.gz$' ${outputs} >> ${hardFilteredVcfs}
      grep 'hard_filtered_with_genotypes.vcf.gz.tbi$' ${outputs} >> ${hardFilteredVcfIndices}
      set +e
      fingerprintingVcf=$(grep 'fingerprinting.vcf.gz$' ${outputs})
      fingerprintingVcfIndex=$(grep 'fingerprinting.vcf.gz.tbi$' ${outputs})
      set -e
      if [[ ! -z "${fingerprintingVcf}" ]]; then
        echo ${fingerprintingVcf} >> ${fingerprintingVcfs}
        echo ${fingerprintingVcfIndex} >> ${fingerprintingVcfIndices}
      fi
      rm ${outputs}
    fi
  done < <(sed -n '1d;p' ${WORKFLOWS_CSV})

  gsutil cp ${sitesOnlyVcfs} ${PART_ONE_FOFN_DIR}/sites_only_vcfs_fofn
  gsutil cp ${sitesOnlyVcfIndices} ${PART_ONE_FOFN_DIR}/sites_only_vcf_indices_fofn
  gsutil cp ${annotationVcfs} ${PART_ONE_FOFN_DIR}/annotation_vcfs_fofn
  gsutil cp ${annotationVcfIndices} ${PART_ONE_FOFN_DIR}/annotation_vcf_indices_fofn
  gsutil cp ${fingerprintingVcfs} ${PART_ONE_FOFN_DIR}/fingerprinting_vcfs_fofn
  gsutil cp ${fingerprintingVcfIndices} ${PART_ONE_FOFN_DIR}/fingerprinting_vcf_indices_fofn
  gsutil cp ${hardFilteredVcfs} ${PART_ONE_FOFN_DIR}/hard_filtered_with_genotypes_vcfs_fofn
  gsutil cp ${hardFilteredVcfIndices} ${PART_ONE_FOFN_DIR}/hard_filtered_with_genotypes_vcf_indices_fofn

  rm ${sitesOnlyVcfs} \
    ${sitesOnlyVcfIndices} \
    ${annotationVcfs} \
    ${annotationVcfIndices} \
    ${fingerprintingVcfs} \
    ${fingerprintingVcfIndices} \
    ${hardFilteredVcfs} \
    ${hardFilteredVcfIndices}
}

# Print the status of a task, color-coded according the status of its shards.
function print_task_status() {
  local -r task=${1} metadataFile=${2} indent=${3}
  tmpStatusesFile=$(mktemp)

  jq -r '.calls["'${task}'"] | group_by(.executionStatus) | map({status: .[0].executionStatus, count: length})[]' ${metadataFile} > ${tmpStatusesFile}

  shardsDone=$(jq -r 'select(.status=="Done") | .count' ${tmpStatusesFile})
  shardsRunning=$(jq -r 'select(.status=="Running") | .count' ${tmpStatusesFile})
  shardsFailed=$(jq -r 'select(.status=="Failed") | .count' ${tmpStatusesFile})
  shardsRetried=$(jq -r 'select(.status=="RetryableFailure") | .count' ${tmpStatusesFile})

  taskStatus=''
  if [[ ${shardsFailed} -eq 0 ]] && [[ ${shardsRunning} -eq 0 ]]; then
    taskStatus='\033[0;32m' # Succeeded
  elif [[ ${shardsRunning} -ne 0 ]] && [[ ${shardsFailed} -ne 0 ]]; then
    taskStatus='\033[0;33m' # Running but will fail
  elif [[ ${shardsRunning} -ne 0 ]]; then
    taskStatus='\033[0;34m' # Running
  else
    taskStatus='\033[0;31m' # Failed
  fi

  echo -e "${taskStatus}${indent}\t${task}\t${shardsRunning:-0} Running, ${shardsDone:-0} Done, ${shardsRetried:-0} Preempted ${shardsFailed:-0} Failed\033[0m"
  rm ${tmpStatusesFile}

}

# Print the status of a workflow (or sub-workflow)
print_workflow_status() {
  local -r metadataFile=${1} expandSubWorkflows=${2} indent=${3:-""}
  local -r -a tasks=($(jq -r '.calls | keys | .[]' ${metadataFile}))

  for task in ${tasks[@]}; do
    if $(jq '.calls["'${task}'"][0] | has("subWorkflowMetadata")' ${metadataFile}) && ${expandSubWorkflows}; then
      local -r subWorkflowName=${task}
      echo -e "\tSubWorkflow ${subWorkflowName}"

      for i in $(seq 0 $(($(jq -r '.calls["'${subWorkflowName}'"] | length ' ${metadataFile}) - 1))); do
        tmpSubWorkflowMetadataFile=$(mktemp)
        jq '.calls["'${subWorkflowName}'"]['${i}'].subWorkflowMetadata' ${metadataFile} > ${tmpSubWorkflowMetadataFile}
        print_workflow_status ${tmpSubWorkflowMetadataFile} ${expandSubWorkflows} "${indent}\t"
        rm ${tmpSubWorkflowMetadataFile}
      done
    else
      print_task_status ${task} ${metadataFile} ${indent}
    fi
  done
}

# Entrypoint for the `do_monitor` per-chromosome logic
function monitor_status() {
  local -r chromosome=${1} workflowId=${2} projectNum=${3} expandSubWorkflows=${4}
  tmpMetadataFile=$(mktemp)
  do_metadata ${chromosome} ${expandSubWorkflows} true > ${tmpMetadataFile}
  local -r workflowStatus=$(jq -r .status ${tmpMetadataFile})
  local -r workflowName=$(jq -r .workflowName ${tmpMetadataFile})
  echo -e "$(tput bold)$(tput smul)${chromosome}$(tput rmul)\t${workflowId}\t${workflowStatus}\t${COMPUTE_PROJECT_PREFIX}${projectNum}$(tput sgr0)"
  if [[ "${workflowStatus}" != "Succeeded" ]] || ${expandSubWorkflows}; then
    print_workflow_status ${tmpMetadataFile} ${expandSubWorkflows}
  else
    echo -e "${workflowName} is done"
  fi
  rm ${tmpMetadataFile}
}

# Monitor workflows
# Print a pretty log to the screen detailing task statuses of a workflow
function do_monitor() {
  local -r expandSubWorkflows=${1:-false} chr=${2}

  if [[ -z "${chr}" ]]; then
    while IFS=, read chromosome workflowId projectNum failedIds; do
      if [[ ! -z ${workflowId} ]]; then
        monitor_status ${chromosome} ${workflowId} ${projectNum} ${expandSubWorkflows}
        echo
      fi
    done < <(sed -n '1d;p' ${WORKFLOWS_CSV})
  else
    local -r id=$(get_workflow_id_from_log ${chr})
    assert_id_exists ${chr} ${id}
    local -r project=$(grep "^${chr}," ${WORKFLOWS_CSV} | cut -d ',' -f 3)
    monitor_status ${chr} ${id} ${project} ${expandSubWorkflows}
  fi
}

# Get the metadata of a workflow
# Optionally expand sub-workflows or only the metadata necessary for monitoring
function do_metadata() {
  local -r chromosome=${1} expandSubWorkflows=${2:-false} beSuperSlim=${3:-false}
  local -r id=$(get_workflow_id_from_log ${chromosome})
  assert_id_exists ${chromosome} ${id}

  local urlParameters=""
  if ${expandSubWorkflows} && ${beSuperSlim}; then
    urlParameters="?${EXPAND_SUB_WORKFLOWS_ARG}&${KEYS_TO_INCLUDE}"
  elif ${expandSubWorkflows}; then
    urlParameters="?${EXPAND_SUB_WORKFLOWS_ARG}"
  elif ${beSuperSlim}; then
    urlParameters="?${KEYS_TO_INCLUDE}"
  fi

  do_get ${CROMWELL_URL}/${id}/metadata${urlParameters}
}

# Get the status of a chromosome's workflow
function do_status() {
  local -r chromosome=${1}
  local -r id=$(get_workflow_id_from_log ${chromosome})
  assert_id_exists ${chromosome} ${id}
  do_get ${CROMWELL_URL}/${id}/status | jq -r .status
}

# Get the outputs of a chromosomes workflow, or get the outputs of all chromosomes
# Rename and pad shard indices to make things easily sortable by downstream folks
function do_outputs() {
  local -r chr=${1}
  if [[ -z "${chr}" ]]; then
    while IFS=, read chromosome workflowId projectNum failedIds; do
      if [[ ! -z ${workflowId} ]]; then
        do_get ${CROMWELL_URL}/${workflowId}/outputs | jq '.'
      fi
    done < <(sed -n '1d;p' ${WORKFLOWS_CSV})
  else
    local -r id=$(get_workflow_id_from_log ${chr})
    assert_id_exists ${chr} ${id}
    do_get ${CROMWELL_URL}/${id}/outputs  | jq  '.outputs'
  fi
}

# Copy the outputs of part one to a safe location that doesn't have a TTL
function copy_part_one_outputs_to_safe_location() {
  local -r chr_arg=${1:-"chr*"}

  while IFS=, read chromosome workflowId projectNum failedIds; do
    if [[ ${chromosome} == ${chr_arg} ]]; then

      outputDir="${OUTPUT_CLOUD_DIR}/${PART_ONE_OUTPUTS_SUBDIR}/${chromosome}/"
      existingFilesFile=$(mktemp)

      if gsutil -q stat "${outputDir}*"; then
        gsutil ls "${outputDir}**/*" | xargs -n1 basename | sort > ${existingFilesFile}
      else
        >&2 echo "No files have been copied for ${chromosome} yet. Copying all files now."
      fi

      outputs=$(do_outputs ${chromosome})
      siteOnlyVcf=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.output_vcf"]')
      sitesOnlyVcfIndex=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.output_vcf_index"]')
      annotationVcf=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.annotation_db_vcf"]')
      annotationVcfIndex=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.annotation_db_vcf_index"]')
      scatteredIntervals=($(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.output_intervals"][]'))
      outputHardFilteredVcfs=($(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.output_hard_filtered_vcfs"][]'))
      outputHardFilteredVcfIndices=($(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.output_hard_filtered_vcf_indices"][]'))
      genomicsDatabases=($(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartOne.genomics_databases"][]'))
      fingerprintingVcf=$(echo ${outputs} | jq -r 'select(.["JointGenotypingByChromosomePartOne.output_fingerprinting_vcf"]!=null) | .["JointGenotypingByChromosomePartOne.output_fingerprinting_vcf"]')
      fingerprintingVcfIndex=$(echo ${outputs} | jq -r 'select(.["JointGenotypingByChromosomePartOne.output_fingerprinting_vcf_index"]!=null) | .["JointGenotypingByChromosomePartOne.output_fingerprinting_vcf_index"]')

      for gdb in ${genomicsDatabases[@]}; do
        (
        fileName=$(basename ${gdb})
        dirName=$(basename $(dirname ${gdb}))
        if [[ ${dirName} == chr2 ]]; then
          dirName="shard-1"
        fi
        paddedName="shard-$(printf "%04d" $(echo ${dirName} | sed -E 's/.*-([0-9]+)/\1/g'))"
        copyName="${paddedName}.${fileName}"

        if grep -q "${copyName}$" ${existingFilesFile}; then
          >&2 echo "Skipping copying of ${gdb} as it already exists in the destination directory"
          continue
        else
          gsutil cp ${gdb} ${outputDir}${copyName}
        fi
        ) &
      done
      wait

      for hardFilteredVcf in ${outputHardFilteredVcfs[@]}; do
        (
        fileName=$(basename ${hardFilteredVcf})
        dirName=$(basename $(dirname ${hardFilteredVcf}))
        shardIndex=$(echo ${fileName} | sed -E 's|.*\.([0-9]+)\.variant_filtered.*|\1|g')
        paddedShardIndex=$(printf "%04d" ${shardIndex})
        paddedName=$(echo ${fileName} | sed -E 's|(.*\.)([0-9]+)(\.variant_filtered.*)|\1'${paddedShardIndex}'\3|g')
        copyName=$(echo ${paddedName} | sed 's|variant_filtered|hard_filtered_with_genotypes|g')

        if grep -q "${copyName}$" ${existingFilesFile}; then
          >&2 echo "Skipping copying of ${hardFilteredVcf} as it already exists in the destination directory"
          continue
        else
          gsutil cp ${hardFilteredVcf} ${outputDir}${copyName}
        fi
        # Check hardFilteredVcf index separately
        if grep -q "${copyName}.tbi$" ${existingFilesFile}; then
          >&2 echo "Skipping copying of ${hardFilteredVcf}.tbi as it already exists in the destination directory"
          continue
        else
          gsutil cp ${hardFilteredVcf}.tbi ${outputDir}${copyName}.tbi
        fi
        ) &
      done
      wait

      allFilesFile=$(mktemp)
      allFiles=(${siteOnlyVcf} ${sitesOnlyVcfIndex} ${annotationVcf} ${annotationVcfIndex} ${scatteredIntervals[@]})
      if [[ ! -z "${fingerprintingVcf}" ]]; then
        allFiles+=(${fingerprintingVcf} ${fingerprintingVcfIndex})
      fi
      printf '%s\n' "${allFiles[@]}" | sort > ${allFilesFile}
      fileNamesToCopy=($(comm -23 <(cat ${allFilesFile} | xargs -n1 basename | sort) ${existingFilesFile}))

      set +e
      filesToCopy=()
      for fileName in ${fileNamesToCopy[@]}; do
        fullFilePath=$(grep "${fileName}$" ${allFilesFile})
        if [[ ! -z "${fullFilePath}" ]]; then
          filesToCopy+=(${fullFilePath})
        fi
      done
      set -e

      if [[ ${#filesToCopy[@]} -gt 0 ]]; then
        >&2 echo "copying ${#filesToCopy[@]} files for ${chromosome}"
        gsutil -m cp ${filesToCopy[@]} ${outputDir}
      fi

      rm ${allFilesFile}
      rm ${existingFilesFile}
    fi
  done < <(sed -n '1d;p' ${WORKFLOWS_CSV})
}

# Copy the outputs of part two to a safe location that doesn't have a TTL
function copy_part_two_outputs_to_safe_location() {
  local -r outputs=$(do_outputs "part-two")
  local -r -a outputVcfs=($(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartTwo.output_vcfs"][]'))
  local -r -a outputVcfIndices=($(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartTwo.output_vcf_indices"][]'))
  local -r annotationVcf=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartTwo.annotation_db_vcf"]')
  local -r annotationVcfIndex=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartTwo.annotation_db_vcf_index"]')
  local -r summaryMetrics=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartTwo.summary_metrics_file"]')
  local -r detailMetrics=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartTwo.detail_metrics_file"]')
  local -r fingerprintMetrics=$(echo ${outputs} | jq -r '.["JointGenotypingByChromosomePartTwo.crosscheck_metrics"]')

  gsutil -m cp ${outputVcfs[@]} ${outputVcfIndices[@]} ${annotationVcf} ${annotationVcfIndex} ${summaryMetrics} ${detailMetrics} ${fingerprintMetrics} "${OUTPUT_CLOUD_DIR}/part_two_outputs/"
}

# Get all the workflows with a given status
function get_workflows_with_status() {
  local -r status=${1}
  while IFS=, read chromosome workflowId projectNum failedIds; do
    if [[ ! -z ${workflowId} ]]; then
      workflowStatus=$(do_status ${chromosome})
      if [[ "${workflowStatus}" == "${status}" ]]; then
        echo ${chromosome}
      fi
    fi
  done < <(sed -n '1d;p' ${WORKFLOWS_CSV})
}

# Get all the workflow IDs for a chromosome
function get_all_workflow_ids() {
  local -r chr=${1}
  if [[ -z "${chr}" ]]; then
    while IFS=, read chromosome workflowId projectNum previousIds; do
      if [[ ! -z ${workflowId} ]]; then
        echo ${workflowId} $(echo ${previousIds} | sed 's/|/ /g')
        echo
      fi
    done < <(sed -n '1d;p' ${WORKFLOWS_CSV})
  else
    local -r workflowId=$(grep "^${chr}," ${WORKFLOWS_CSV} | cut -d ',' -f 2)
    local -r -a previousIds=($(grep "^${chr}," ${WORKFLOWS_CSV} | cut -d ',' -f 4 | sed 's/|/ /g'))
    local -a reversedPreviousIds=()
    for (( idx=${#previousIds[@]}-1 ; idx>=0 ; idx-- )) ; do
      reversedPreviousIds+=(${previousIds[idx]})
    done
    echo ${workflowId} ${reversedPreviousIds[@]}
  fi
}

# Calculate the cost of running a chromosome
# Takes into account all previous runs
function get_chromosome_cost() {
  local -r chr=${1}
  local -r -a workflowIds=($(get_all_workflow_ids ${chr}))
  local -a totalCosts=()
  for workflowId in ${workflowIds[@]}; do
    tmpMetadata=$(mktemp)
    do_get "${CROMWELL_URL}/${workflowId}/metadata?${EXPAND_SUB_WORKFLOWS_ARG}&${COST_KEYS_TO_EXCLUDE}" > ${tmpMetadata}
    costScriptOutput=$(${DSDE_PIPELINES_DIR}/scripts/calculate_cost.py -m ${tmpMetadata} --only_total || echo 0)
    >&2 echo "${chr} ${workflowId}: ${costScriptOutput}"
    rawCost=$(echo ${costScriptOutput} | cut -d ' ' -f 3 | sed 's/[eE]+*/\*10\^/')
    totalCosts+=(${rawCost})
    rm ${tmpMetadata}
    # Let cromwell recover...
    sleep 30
  done
  sum=$( IFS="+"; bc <<< "${totalCosts[*]}" )
  echo "${chr} Total Cost: ${sum}"
}

# Get the total cost of this a project
# This can cause a strain on Cromwell, so be careful
function get_total_cost() {
  local -a chromosomeCosts=()
  >&2 echo "This will take a long time..."
  while IFS=, read chromosome workflowId projectNum previousIds; do
      chromosomeCost=$(get_chromosome_cost ${chromosome})
      >&2 echo ${chromosomeCost}
      rawChromosomeCost=$(echo ${chromosomeCost} | cut -d ' ' -f 4)
      chromosomeCosts+=(${rawChromosomeCost})
  done < <(sed -n '1d;p' ${WORKFLOWS_CSV})
  sum=$( IFS="+"; bc <<< "${chromosomeCosts[*]}" )
  echo "Joint Genotyping Total Cost: ${sum}"
}

# Clear a chromosome's workflow so another one can be run
# Don't abort the workflow, just let it run
# Also useful to clear a previously failed workflow from the log
function void_workflow() {
  local -r chromosome=${1}

  local -r id=$(get_workflow_id_from_log ${chromosome})
  assert_id_exists ${chromosome} ${id}
  sed -i.bak -E 's/^'${chromosome}',([a-f0-9-]+),([0-9]+),([|a-f0-9-]*)/'${chromosome}',,\2,\3\1|/g' ${WORKFLOWS_CSV}
  rm ${WORKFLOWS_CSV}.bak
}

# Abort the currently running workflow of a chromosome.
function do_abort() {
  local -r chromosome=${1}
  local -r id=$(get_workflow_id_from_log ${chromosome})
  assert_id_exists ${chromosome} ${id}
  do_post ${CROMWELL_URL}/${id}/abort | jq .
  void_workflow ${chromosome}
}

expandSubWorkflows="false"
chr=""
projectNum=$(( ( RANDOM % ${PROJECT_NUMBER_RANGE[1]} )  + 1 ))
command=${1}
shift

if [[ ${command} == "-h" ]]; then
  >&2 echo '`submit-part-one -c chrN [-p P]` - Submit `chrN` to cromwell for part one of joint genotyping. Optionally, a number `P` can be given to specify a google project to run on. If not given, the script will choose a random project.'
  >&2 echo '`submit-part-two [-p P]` - Submit the results of part-one to cromwell for part two of joint genotyping. Optionally, a number `P` can be given to specify a google project to run on. If not given, the script will choose a random project.'
  >&2 echo '`status -c chrN` - Get the status of the current workflow for `chrN`.'
  >&2 echo '`monitor [-x] [-c chrN]` - Monitor currently running workflows for all chromosomes. Providing `-x` will expand sub-workflow monitoring. Providing `-c chrN` will monitor only the provided chromosome instead of all running ones.'
  >&2 echo '`metadata -c chrN [-x]` - Get the metadata of the current workflow for `chrN`. Providing `-x` will expand sub-workflow metadata.'
  >&2 echo '`abort -c chrN` - Abort the workflow for `chrN` so another workflow for `chrN` can be run.'
  >&2 echo '`void -c chrN` - Allow another workflow to be run for `chrN`, but do not abort the current workflow.'
  >&2 echo '`outputs -c chrN` - Get the outputs of a successful `chrN` workflow.'
  >&2 echo '`copy-part-one-outputs [-c chrN]` - Copy the outputs of part one out of the execution directory and into a safe location. Providing `-c chrN` will only copy the specified chromosomes outputs.'
  >&2 echo '`copy-part-two-outputs` - Copy the outputs of part two out of the execution directory and into a safe location.'
  >&2 echo '`get-failed` - Get all the chromosomes that have failed.'
  >&2 echo '`get-succeeded` - Get all the chromosomes that have succeeded.'
  >&2 echo '`get-running` - Get all the chromosomes that are still running.'
  >&2 echo '`chromosome-cost -c chrN` - Calculate the cost of running `chrN`.'
  >&2 echo '`total-cost` - Calculate the total cost of all workflows run.'
  exit
fi

while getopts "c:p:x" option; do
    case ${option} in
      c)  chr=${OPTARG} ;;
      x)  expandSubWorkflows="true" ;;
      p)
        projectNum=${OPTARG}
        if [[ ${projectNum} -gt ${PROJECT_NUMBER_RANGE[0]} ]] || [[ ${projectNum} -lt ${PROJECT_NUMBER_RANGE[0]} ]]; then
          >&2 echo "Invalid project number given. Please provide a project in the range ${PROJECT_NUMBER_RANGE[0]} - ${PROJECT_NUMBER_RANGE[0]} inclusive."
          exit 1
        fi
      ;;
  esac
done
shift $((OPTIND -1))

if [[ -z ${chr} ]]
then
  declare -r -a non_chr_commands=(monitor get-failed get-succeeded get-running total-cost submit-part-two copy-part-one-outputs copy-part-two-outputs setup)
  case "${non_chr_commands[@]}" in
    *"${command}"*) true;;
    *)
      >&2 echo "${command} requires an argument declared with [-c chrN]" 2>&1
      exit 1
      ;;
  esac
fi

case ${command} in
  submit-part-one) do_submit_part_one ${chr} ${projectNum} ;;
  submit-part-two) do_submit_part_two ${projectNum} ;;
  status) do_status ${chr} ;;
  abort) do_abort ${chr} ;;
  void) void_workflow ${chr} ;;
  metadata) do_metadata ${chr} ${expandSubWorkflows} ;;
  monitor) do_monitor ${expandSubWorkflows} ${chr} ;;
  outputs) do_outputs ${chr} ;;
  copy-part-one-outputs) copy_part_one_outputs_to_safe_location ${chr} ;;
  copy-part-two-outputs) copy_part_two_outputs_to_safe_location ;;
  get-failed) get_workflows_with_status "Failed" ;;
  get-succeeded) get_workflows_with_status "Succeeded" ;;
  get-running) get_workflows_with_status "Running" ;;
  chromosome-cost) get_chromosome_cost ${chr} ;;
  total-cost) get_total_cost ;;
  *)
    >&2 echo "${command} is not a valid command"
    exit 1
  ;;
esac
