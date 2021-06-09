#!/usr/bin/env bash
set -e

VERIFY_BAM_ID_COMMIT="c1cba76e979904eb69c31520a0d7f5be63c72253"
TIMESTAMP=$(date -u +"%Y-%m-%d_%H-%M-%SZ")

TAG="$VERIFY_BAM_ID_COMMIT-$TIMESTAMP"
GCR=us.gcr.io/broad-gotc-prod/verify-bam-id

docker build -t ${GCR}:${TAG} \
		--build-arg GIT_HASH="$VERIFY_BAM_ID_COMMIT" .

docker -- push ${GCR}:${TAG}

echo - e "$GCR:$TAG" >> ./build_docker_version.tsv