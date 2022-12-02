#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.7
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/imputation-bcf-vcf"
# QUAY_URL="quay.io/broadinstitute/gotc-prod-imputation_bcf_vcf"

#BCFTOOLS version
BCFTOOLS_VERSION="1.10.2"

#VCFTOOLS version
VCFTOOLS_VERSION="0.1.16"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-b|--bcf] [-v|--vcf] [-t|--tools] -- script to build the Imputation Bcf/Vcf tools image and push to GCR & Quay

where:
    -h|--help Show help text
    -b|--bcf Version of BCFTOOLS to use (default: BCFTOOLS_VERSION=${BCFTOOLS_VERSION})
    -v|--vcf Version of VCFTOOLS to use (default: VCFTOOLS_VERSION=${VCFTOOLS_VERSION})
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
        -b|--bcf)
        BCFTOOLS_VERSION="$2"
        shift
        shift
        ;;
        -v|--vcf)
        VCFTOOLS_VERSION="$2"
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$BCFTOOLS_VERSION-$VCFTOOLS_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg BCFTOOLS_VERSION="$BCFTOOLS_VERSION" \
        --build-arg VCFTOOLS_VERSION="$VCFTOOLS_VERSION" \
        --no-cache $DIR   
    docker push "$GCR_URL:$IMAGE_TAG"

#    echo "tagging and pushing Quay Image"
#    docker tag "$GCR_URL:$IMAGE_TAG" "$QUAY_URL:$IMAGE_TAG"
#    docker push "$QUAY_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"