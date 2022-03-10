#!/usr/bin/env bash
set -e

# VerifyBamID is based off a specific githash and unlikely to change
VERIFY_BAM_ID_VERSION="c1cba76e979904eb69c31520a0d7f5be63c72253"
DOCKER_IMAGE_VERSION="1.0.1"
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/verify-bam-id"
QUAY_URL="quay.io/broadinstitute/gotc-prod-verify_bam_id"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-v|--version] [-t|tools] -- script to build the VerifyBamID image and push to GCR & Quay

where:
    -h|--help Show help text
    -v|--version Git hash of the VerifyBamID version to use (default: $VERIFY_BAM_ID_VERSION)
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
        VERIFY_BAM_ID_VERSION="$2"
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$VERIFY_BAM_ID_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build --no-cache -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg GIT_HASH="$VERIFY_BAM_ID_VERSION" "$DIR"
    docker push "$GCR_URL:$IMAGE_TAG"

    echo "tagging and pushing Quay Image"
    docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
    docker push "$QUAY_URL:$IMAGE_TAG"

    echo "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"