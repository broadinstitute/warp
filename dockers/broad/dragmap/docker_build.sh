#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.1.2
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/dragmap"
QUAY_URL="quay.io/broadinstitute/gotc-prod-dragmap"

# DRAGMAP version
DRAGMAP_VERSION="1.2.1"

# PICARD version
PICARD_VERSION="2.26.10"

# SAMTOOLS version
SAMTOOLS_VERSION="1.11"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-d|--dragmap_version] [-p|--picard_version] [-s|--samtools_version] [-t|--tools] -- script to build the DRAGMAP image and push to GCR & Quay

where:
    -h|--help Show help text
    -d|--dragmap_version Version of DRAGMAP to use (default: DRAGMAP=${DRAGMAP_VERSION})
    -p|--picard_version Version of PICARD to use (default: PICARD=${PICARD_VERSION})
    -s|--samtools_version Version of SAMTOOLS to use (default: SAMTOOLS=${SAMTOOLS_VERSION})
    -t|--tools Show tools needed to run script
    "

function main(){
    for t in "${TOOLS[@]}"; do which $t >/dev/null || ok=no; done
        if [[ $ok == no ]]; then
            echo "Missing one of the following tools: "
            for t in "${TOOLS[@]}"; do echo "$t"; done
            exit 1
        fi

    while [[ $# -gt 0 ]]
    do 
    key="$1"
    case $key in
        -d|--dragmap_version)
        DRAGMAP_VERSION="$2"
        shift
        shift
        ;;
        -p|--picard_version)
        PICARD_VERSION="$2"
        shift
        shift
        ;;
        -s|--samtools_version)
        SAMTOOLS_VERSION="$2"
        shift
        shift
        ;;
        -h|--help)
        echo "${HELP}"
        exit 0
        ;;
        -t|--tools)
        for t in "${TOOLS[@]}"; do echo $t; done
        exit 0
        ;;
        *)
        shift
        ;;
    esac
    done

    IMAGE_TAG="${DOCKER_IMAGE_VERSION}-${DRAGMAP_VERSION}-${PICARD_VERSION}-${SAMTOOLS_VERSION}-${TIMESTAMP}"

    echo "building and pushing GCR Image - ${GCR_URL}:${IMAGE_TAG}"
    docker build -t "${GCR_URL}:${IMAGE_TAG}" \
        --build-arg DRAGMAP_VERSION="${DRAGMAP_VERSION}" \
        --build-arg PICARD_VERSION="${PICARD_VERSION}" \
        --build-arg SAMTOOLS_VERSION="${SAMTOOLS_VERSION}" \
        --no-cache $DIR
    docker push "${GCR_URL}:${IMAGE_TAG}"

    echo "tagging and pushing Quay Image"
    docker tag "${GCR_URL}:${IMAGE_TAG}" "${QUAY_URL}:${IMAGE_TAG}"
    docker push "${QUAY_URL}:${IMAGE_TAG}"

    echo -e "${GCR_URL}:${IMAGE_TAG}" >> "${DIR}/docker_versions.tsv"
    echo "done"
}

main "$@"
