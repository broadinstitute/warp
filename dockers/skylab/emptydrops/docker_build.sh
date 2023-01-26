#!/bin/bash
#fail-fast
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.1
TIMESTAMP=$(date +"%s")
DIR=$(cd "$(dirname "$0")" && pwd)
TAG=$1
CACHING=$2

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/empty-drops"

#R Version
R_VERSION="4.2.2"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-t|tools] [TAG] [CACHING] -- script to build the empty-drops image and push to GCR
where:
    -h|--help Show help text
    -s|--r_version Version of R to use (default: $R_VERSION)
    -t|--tools Show tools needed to run script
    TAG : Set a fixed tag string (this will override our tagging convention). Useful for development or on CI/CD to reduce wasted images
    CACHING : Set to ON to enable caching. Any other value or simply not setting it will default to no-cache builds.
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
        -s|--R_VERSION)
        R_VERSION="$2"
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
    if [ -z $tag ]; 
    then
        IMAGE_TAG="$TAG"
    else
        IMAGE_TAG="$DOCKER_IMAGE_VERSION-$R_VERSION-$TIMESTAMP"
    fi

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"

    if [ $CACHING == "ON" ]; then
        docker build -t "$GCR_URL:$IMAGE_TAG" --build-arg R_VERSION="$R_VERSION" "$DIR"
    else
        docker build -t "$GCR_URL:$IMAGE_TAG" --build-arg R_VERSION="$R_VERSION" --no-cache "$DIR"
    fi
    docker push "$GCR_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"