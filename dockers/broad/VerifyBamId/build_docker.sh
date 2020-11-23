#!/usr/bin/env bash

DIR="$( cd "$( dirname "$0" )" && pwd )"

VERIFY_BAM_ID_COMMIT="c1cba76e979904eb69c31520a0d7f5be63c72253"
TIMESTAMP=$(date +"%s")
DOCKER_TAG="$VERIFY_BAM_ID_COMMIT-$TIMESTAMP"

docker build -t us.gcr.io/broad-gotc-prod/verify-bam-id:${DOCKER_TAG} --build-arg GIT_HASH="$VERIFY_BAM_ID_COMMIT" "$DIR"
gcloud docker -- push us.gcr.io/broad-gotc-prod/verify-bam-id:${DOCKER_TAG}
