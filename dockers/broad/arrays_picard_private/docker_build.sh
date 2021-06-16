#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=4.0.10
TIMESTAMP=$(date -u +"%Y-%m-%d")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL=us.gcr.io/broad-arrays-prod/arrays-picard-private
IMAGE_TAG="$DOCKER_IMAGE_VERSION-$TIMESTAMP"

# Picard private artifact
PICARD_PRIVATE_VERSION=61af9bff4587783e5981a496f422ea36102482b5
ARTIFACTORY_URL="https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/picard-private/"

# Necessary tools and help text
TOOLS=(docker gcloud vault jq)
HELP="$(basename "$0") [-h|--help] [-v|--version] [-t|--tools] -- script to build the picard private image and push to GCR

where:
    -h|--help Show help text
    -v|--version Git hash of the picard private version to use (default: $PICARD_PRIVATE_VERSION)
    -t|--tools Show tools needed to run script
    "


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
    PICARD_PRIVATE_VERSION="$2"
    shift
    shift
    ;;
    -h|--help)
    echo "$HELP"
    exit 0
    ;;
    -t|--tools)
    for t in "${TOOLS[@]}"; do echo $t; done
    shift
    ;;
    *)
    shift
    ;;
esac
done

echo "downloading picard private jar -  $ARTIFACTORY_URL/$PICARD_PRIVATE_VERSION"
curl "$ARTIFACTORY_URL/$PICARD_PRIVATE_VERSION" > "$DIR/picard-private.jar"

echo "building & pushing GCR Image - $GCR_URL:$IMAGE_TAG"
docker build --no-cache -t "$GCR_URL:$IMAGE_TAG" .
docker push $GCR_URL:$IMAGE_TAG

echo "removing picard private jar - $DIR/picard-private.jar"
rm "$DIR/picard-private.jar"

echo "$GCR_URL:$IMAGE_TAG $PICARD_PRIVATE_VERSION" >> "$DIR/docker_versions.tsv"

echo "done"



