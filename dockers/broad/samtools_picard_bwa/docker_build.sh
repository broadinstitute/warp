#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.2
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/samtools-picard-bwa"
QUAY_URL="quay.io/broadinstitute/gotc-prod-samtools_picard_bwa"

# BWA version
BWA_VERSION="0.7.15"

# PICARD PUBLIC version
PICARD_PUBLIC_VERSION="2.26.10"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-b|--bwa_version] [-p|--picard_public_version] [-t|tools] -- script to build the SAMTOOLS/PICARD/BWA image and push to GCR & Quay

where:
    -h|--help Show help text
    -b|--bwa_version Version of BWA to use (default: BWA=$BWA_VERSION)
    -p|--picard_public_version Version of PICARD_PUBLIC to use (default: PICARD_PUBLIC=$PICARD_PUBLIC_VERSION)
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
        -b|--bwa_version)
        BWA_VERSION="$2"
        shift
        shift
        ;;
        -p|--picard_public_version)
        PICARD_PUBLIC_VERSION="$2"
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$BWA_VERSION-$PICARD_PUBLIC_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg BWA_VERSION="$BWA_VERSION" \
        --build-arg PICARD_PUBLIC_VERSION="$PICARD_PUBLIC_VERSION" \
        --no-cache $DIR
    docker push "$GCR_URL:$IMAGE_TAG"

    echo "tagging and pushing Quay Image"
    docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
    docker push "$QUAY_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"
