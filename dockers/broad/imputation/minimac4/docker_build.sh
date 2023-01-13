#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.6
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/imputation-minimac4"
#QUAY_URL="quay.io/broadinstitute/gotc-prod-imputation_minimac"

# Minimac4 version
MINIMAC4_VERSION="1.0.2"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-m|--minimac] [-t|--tools] -- script to build the Imputation Minimac4 image and push to GCR & Quay

where:
    -h|--help Show help text
    -m|--minimac Version of Minimac4 to use (default: MINIMAC4_VERSION=${MINIMAC4_VERSION})
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
        -m|--minimac)
        MINIMAC_VERSION="$2"
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$MINIMAC4_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg MINIMAC4_VERSION="$MINIMAC4_VERSION" \
        --no-cache $DIR   
    docker push "$GCR_URL:$IMAGE_TAG"

#    echo "tagging and pushing Quay Image"
#    docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
#    docker push "$QUAY_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"