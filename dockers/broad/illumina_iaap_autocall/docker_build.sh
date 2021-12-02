#!/bin/bash
set -e

# Update verson when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.2
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags"
GCR_URL="us.gcr.io/broad-gotc-prod/illumina-iaap-autocall"
QUAY_URL="quay.io/broadinstitute/gotc-prod-illumina_iaap_autocall"


# Iaap cli version
IAAP_CLI_VERSION="iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-t|--tools] -- script to build the Illumina IAAP image and push to GCR & Quay

where:
    -h|--help Show help text
    -v|--version Zip version of Zcall to use (default: $IAAP_CLI_VERSION)
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
        IAAP_CLI_VERSION="$2"
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

    VERSION_NUMBER=$(echo $IAAP_CLI_VERSION | grep -Eo '[-][0-9]+([.][0-9]+)+[-]' | sed 's/-//g' )
    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$VERSION_NUMBER-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build --no-cache -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg IAAP_CLI_VERSION="$IAAP_CLI_VERSION" "$DIR"
    docker push "$GCR_URL:$IMAGE_TAG"

    echo "tagging and pushing Quay Image"
    docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
    docker push "$QUAY_URL:$IMAGE_TAG"

    echo "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"