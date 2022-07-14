#!/bin/bash

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
SNAPTOOLS_VERSION="1.2.3"

GCR_URL="us.gcr.io/broad-gotc-prod/snaptools"
#QUAY_URL="quay.io/broadinstitute/gotc-prod-samtools"


IMAGE_TAG="$DOCKER_IMAGE_VERSION-$SNAPTOOLS_VERSION-$TIMESTAMP"

echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
docker build --no-cache -t "$GCR_URL:$IMAGE_TAG" .
docker push "$GCR_URL:$IMAGE_TAG"

#echo "tagging and pushing Quay Image"
#docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
#docker push "$QUAY_URL:$IMAGE_TAG"

echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"

#docker build . -t quay.io/humancellatlas/snaptools:0.0.1

#echo You can now push to quay.io using the following command
#echo    docker push quay.io/humancellatlas/snaptools:0.0.1