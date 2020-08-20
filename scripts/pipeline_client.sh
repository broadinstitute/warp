#!/usr/bin/env bash

ROOT_DIR=$(git rev-parse --show-toplevel)

CROMWELL_ENV="dev"
EXPAND_SUB_WORKFLOWS="false"

RELEASE_DIR=''

function do_release() {
  local -r wdl=${1}
  RELEASE_DIR=$(mktemp -d)

  ${ROOT_DIR}/scripts/build_workflow_release.sh \
    ${wdl} \
    ${RELEASE_DIR} \
    $(git rev-parse --short HEAD) \
    ${CROMWELL_ENV} > /dev/null

  echo Released to temporary directory ${RELEASE_DIR}
}

function submit_job() {
  local -r wdlName=${1} inputs=${2}

  cromshell submit \
    ${RELEASE_DIR}/${wdlName}/${wdlName}.wdl \
    ${inputs} \
    ${RELEASE_DIR}/${wdlName}/${wdlName}.options.json \
    ${RELEASE_DIR}/dependencies/workflow_dependencies.zip

  rm -rf ${RELEASE_DIR}
}

function do_submit() {
  local -r wdl=${1} inputs=${2}
  do_release ${wdl}
  local -r wdlName=$(basename ${wdl} .wdl)
  submit_job ${wdlName} ${inputs}
}

function do_metadata() {
  local -r workflowId=${1}
  cromshell metadata ${workflowId}
}

function do_status() {
  local -r workflowId=${1}
  cromshell status ${workflowId}
}

function do_abort() {
  local -r workflowId=${1}
  cromshell abort ${workflowId}
}

function usage() {
  >&2 echo "\`pipeline_client.sh\` is a wrapper script around cromshell."
  >&2 echo "It will submit WDLs in dsde-pipelines while also managing dependencies"
  >&2 echo "Usage:"
  >&2 echo "pipeline_client.sh submit \\"
  >&2 echo "  pipelines/dna_seq/germline/wgs/WholeGenomeGermlineSingleSample.wdl \\"
  >&2 echo "  pipelines/dna_seq/germline/wgs/input_files/WholeGenomeGermlineSingleSample.inputs.plumbing.json"
  >&2 echo
  >&2 echo "By default, \`pipeline_client.sh submits to the GoTC Dev cromwell.\`"
  >&2 echo "The cromwell server can be set by supplying \`-e ENV\` after a command."
  >&2 echo "ENV can by one of [dev|staging|prod|jgdev|jgprod]"
  >&2 echo "Example:"
  >&2 echo "pipeline_client.sh submit -e prod \\"
  >&2 echo "  pipelines/dna_seq/germline/wgs/WholeGenomeGermlineSingleSample.wdl \\"
  >&2 echo "  pipelines/dna_seq/germline/wgs/input_files/WholeGenomeGermlineSingleSample.inputs.plumbing.json"
  >&2 echo
  >&2 echo "Available commands are:"
  >&2 echo "  submit"
  >&2 echo "  metadata"
  >&2 echo "  status"
  >&2 echo "  abort"
  >&2 echo "  help"
  >&2 echo
}

COMMAND=${1}
shift

while getopts "e:" option; do
    case ${option} in
      e)  CROMWELL_ENV=${OPTARG} ;;
      *)
        usage
        >&2 echo "Illegal option"
        exit 1
        ;;
  esac
done
shift $((OPTIND -1))

case ${CROMWELL_ENV} in
  dev) export CROMWELL_URL="https://cromwell.gotc-dev.broadinstitute.org" ;;
  staging) export CROMWELL_URL="https://cromwell.gotc-staging.broadinstitute.org" ;;
  prod) export CROMWELL_URL="https://cromwell.gotc-prod.broadinstitute.org" ;;
  jgdev) export CROMWELL_URL="https://cromwell-jg.gotc-dev.broadinstitute.org" ;;
  jgprod) export CROMWELL_URL="https://cromwell-jg.gotc-prod.broadinstitute.org" ;;
  *)
    >&2 echo "Cromwell environment must be one of [dev|staging|prod|jgdev|jgprod]"
    exit 1
    ;;
esac

case ${COMMAND} in
  submit) do_submit ${1} ${2} ;;
  metadata) do_metadata ${1} ;;
  status) do_status ${1} ;;
  abort) do_abort ${1} ;;
  help) usage ;;
  *)
    usage
    >&2 echo "${COMMAND} is not a valid command"
    exit 1
  ;;
esac
