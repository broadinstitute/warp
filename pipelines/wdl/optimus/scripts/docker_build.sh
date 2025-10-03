#!/bin/bash

set -e

DOCKER_IMAGE_TAG=3.0
DIR=$(cd $(dirname $0) && pwd)

QUAY_URL="quay.io/humancellatlas/hca_preprocessing"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-i|image] [-t|tools] -- script to build the HCA preprocessing image
where:
    -h|--help Show help text
    -i|--image Tag for the docker image
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
        -i|--image)
        DOCKER_IMAGE_TAG="$2"
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


    echo "Building and pushing Quay Image"
    docker build --no-cache -t "$QUAY_URL:$DOCKER_IMAGE_TAG" $DIR
    docker push "$QUAY_URL:$DOCKER_IMAGE_TAG"

    echo -e "$QUAY_URL:$DOCKER_IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"