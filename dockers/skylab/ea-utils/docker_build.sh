#!/bin/bash

set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/ea-utils"
#QUAY_URL="quay.io/humancellatlas/ea-utils"

# ea-utils version
EAUTILS_VERSION="1.04.807"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-t|--tools] -- script to build the ea-utils image and push to GCR & Quay

where:
    -h|--help Show help text
    -v|--version Version of ea-utils to use (default: $EAUTILS_VERSION)
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
        -v|--version)
        EAUTILS_VERSION="$2"
        shift
        shift
        ;;
        -h|--help)
        echo "$HELP"
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$EAUTILS_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build --no-cache -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg EAUTILS_VERSION="$EAUTILS_VERSION" "$DIR"
    docker push "$GCR_URL:$IMAGE_TAG"

    #echo "tagging and pushing Quay Image"
    #docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
    #docker push "$QUAY_URL:$IMAGE_TAG"

    echo "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"
