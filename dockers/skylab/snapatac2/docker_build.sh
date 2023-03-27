#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.2
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/"
QUAY_URL="quay.io/broadinstitute/gotc-prod-snapatac2/"

# SnapATAC2 version
SNAPATAC2_VERSION="2.2.0"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-b|--snap_version] [-t|tools] -- script to build the SAMTOOLS/BWA image and push to GCR & Quay

where:
    -h|--help Show help text
    -b|--snap_version Version of SNAPATAC2 to use (default: SNAPATAC2=$SNAPATAC2_VERSION)
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
        -b|--snap_version)
        SNAPATAC2_VERSION="$2"
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

    IMAGE_TAG="snapatac2:$DOCKER_IMAGE_VERSION-$SNAPATAC2_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL$IMAGE_TAG"
    docker build -t "$GCR_URL$IMAGE_TAG" \
        --build-arg SNAPATAC2_VERSION="$SNAPATAC2_VERSION" \
        --no-cache $DIR
    docker push "$GCR_URL$IMAGE_TAG"

    #echo "tagging and pushing Quay Image"
    #docker tag "$GCR_URL$IMAGE_TAG" "$QUAY_URL$IMAGE_TAG"
    #docker push "$QUAY_URL$IMAGE_TAG"

    echo -e "$GCR_URL$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"
