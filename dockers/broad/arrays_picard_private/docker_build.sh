#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=4.0.10
TIMESTAMP=$(date -u +"%Y-%m-%d")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL=us.gcr.io/broad-arrays-prod/arrays-picard-private
IMAGE_TAG=$DOCKER_IMAGE_VERSION-$TIMESTAMP

# Picard private artifact
PICARD_PRIVATE_VERSION=61af9bff4587783e5981a496f422ea36102482b5
ARTIFACTORY_URL=https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/picard-private/$PICARD_PRIVATE

# Necessary tools and help text
TOOLS=(docker gcloud vault jq)
HELP="$(basename "$0") [-h|--help] [-v|--version] [-t|--tools] -- script to build the picard-private image and push to GCR & Dockerhub

where:
    -h Show help text
    -v|--version Git hash of the picard private version to use (default: $PICARD_PRIVATE_VERSION)
    -t|--tools Show tools needed to run script
    "


for t in "${TOOLS[@]}"; do which $t >/dev/null || ok=no; done
    if [ "$ok" == "no" ]; then
        echo "Missing one of the following tools: "
        for t in "${TOOLS[@]}"; do echo "$t" 1>&2; done
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
    exit 1
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

echo "Downloading picard private jar: $ARTIFACTORY_URL"
curl $ARTIFACTORY_URL > $DIR/picard-private.jar

echo "Building & pushing GCR Image: $GCR_URL:$IMAGE_TAG"
docker build --no-cache -t $GCR_URL:$IMAGE_TAG .
#docker push $GCR_URL:$IMAGE_TAG

echo "Removing picard private jar"
rm $DIR/picard-private.jar

echo "$GCR_URL:$IMAGE_TAG\t$PICARD_PRIVATE_VERSION" >> $DIR/docker_versions.tsv
echo "Done"