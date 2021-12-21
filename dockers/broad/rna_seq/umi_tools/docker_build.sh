#!/bin/bash
set -e

# Update verson when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/umi_tools"
QUAY_URL="quay.io/broadinstitute/gotc-prod-umi_tools"

# umi_tools version
UMI_TOOLS_VERSION="1.1.1"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-t|--tools] -- script to build the UMI_TOOLS image and push to GCR & Quay

where:
    -h|--help Show help text
    -v|--version of umi_tools to use (default: $FGBIO_VERSION)
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
        UMI_TOOLS_VERSION="$2"
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$UMI_TOOLS_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build --no-cache -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg UMI_TOOLS_VERSION="$UMI_TOOLS_VERSION" "$DIR"
    docker push "$GCR_URL:$IMAGE_TAG"

    #echo "tagging and pushing Quay Image"
    #docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
    #docker push "$QUAY_URL:$IMAGE_TAG"

    echo "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"