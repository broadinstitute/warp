#!/bin/bash
#fail-fast
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
DIR=$(cd "$(dirname "$0")" && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/build-indices"
#QUAY_URL="quay.io/humancellatlas/build-indices"
#Image path before getting renamed: quay.io/humancellatlas/snss2-indices:1.2.0

#STAR version
STAR_VERSION="2.7.10a"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-t|tools] -- script to build the build-indices image and push to GCR
where:
    -h|--help Show help text
    -s|--star_version Version of STAR to use (default: $STAR_VERSION)
    -t|--tools Show tools needed to run script
    "

function main(){
    for t in "${TOOLS[@]}"; do which "$t" >/dev/null || ok=no; done
        if [[ $ok == no ]]; then
            echo "Missing one of the following tools: "
            for t in "${TOOLS[@]}"; do echo "$t"; done
            exit 1
        fi

    while [[ $# -gt 0 ]]
    do
    key="$1"
    case $key in
        -h|--help)
        echo "$HELP"
        exit 0
        ;;
        -s|--star_version)
        STAR_VERSION="$2"
        shift
        shift
        ;;
        -t|--tools)
        for t in "${TOOLS[@]}"; do echo "$t"; done
        exit 0
        ;;
        *)
        shift
        ;;
    esac
    done

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$STAR_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg STAR_VERSION="$STAR_VERSION" \
        --no-cache "$DIR"
    docker push "$GCR_URL:$IMAGE_TAG"

#    echo "tagging and pushing Quay Image"
#    docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
#    docker push "$QUAY_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"