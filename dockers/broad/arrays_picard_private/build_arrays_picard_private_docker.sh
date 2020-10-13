#!/bin/bash
set -e

# Update DOCKER_IMAGE_VERSION after any substantial changes to the
# image this builds.
#
declare -r DOCKER_IMAGE_VERSION=4.0.10

# Update this when there is a new release of picard-private to use as the
# default jar.
#
declare -r PICARD_PRIVATE_VERSION=61af9bff4587783e5981a496f422ea36102482b5

declare -r ARTIFACTORY=https://broadinstitute.jfrog.io/broadinstitute
declare -r LIBS_SNAPSHOT_LOCAL=$ARTIFACTORY/libs-snapshot-local
declare -r LIBS_RELEASE_LOCAL=$ARTIFACTORY/libs-release-local

declare -r PICARD_PRIVATE=$LIBS_SNAPSHOT_LOCAL/org/broadinstitute/picard-private/$PICARD_PRIVATE_VERSION
declare -r PICARD_PRIVATE_URL=$PICARD_PRIVATE/jars/picard-private-all-$PICARD_PRIVATE_VERSION.jar

# Show usage and fail if there is a problem with the command line.
#
function help {
    local scriptName="$1" jar="$3" ok=yes
    local -r -a tools=(docker gcloud vault jq)
    local -r -a usage=(
        ''
        "$scriptName: Build a GotC Arrays docker image and push it to the GCR."
        ''
        'The GCR is the Google Cloud Registry.'
        ''
        "Usage: $scriptName [<jar>]"
        ''
        'Where: <jar> is the path to a picard-private.jar.'
        ''
        "The default <jar> is $PICARD_PRIVATE_URL"
        ''
        "Example: $scriptName ./dist/picard-private.jar"
        ''
        "Note: You need at least these tools installed: ${tools[*]}"
        ''
    )
    for t in "${tools[@]}"; do which $t >/dev/null || ok=no; done
    if test $ok = no
    then
        for line in "${usage[@]}"; do echo "$line" 1>&2; done
        exit 1
    fi
}

# Copy the (possibly remote) picard-private.jar ($jar) binary
# to the local docker build directory ($dir).
#
function copyPicardPrivateToLocalDockerDirectory () {
    local -r jar="$1" dir="$2"
    local picardPrivate="$jar"
    if [ "$picardPrivate" ]; then
        echo scp "$picardPrivate" "$dir"/picard-private.jar
        scp "$picardPrivate" "$dir"/picard-private.jar
    else
        echo "curl $PICARD_PRIVATE_URL > $dir/picard-private.jar"
        curl $PICARD_PRIVATE_URL > "$dir"/picard-private.jar
    fi
}

# Make the Dockerfile.
#
function makeDockerfile () {
    # This base docker image is in turn based on a docker image from Alpine, with glibc added (needed for snappy)
cat << EOF > ./Dockerfile
FROM frolvlad/alpine-glibc

ENV TERM=xterm-256color

WORKDIR /usr/gitc

# Install dependencies.
RUN apk --update add bash curl findutils jq wget openjdk8-jre python2 unzip

# download the gsutil install script
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash -O - | bash

# Set up Vault
ENV NO_VAULT=true
RUN curl https://releases.hashicorp.com/vault/1.0.2/vault_1.0.2_linux_amd64.zip > temp.zip \
    && unzip temp.zip \
    && rm temp.zip \
    && mv vault /usr/local/bin/

COPY . .
EOF
}

# Run docker build, then tag and push the image to the GCR.
#
# By default, set `--no-cache=true` here so we will always build our
# images from scratch instead of possibly missing failing steps in the
# process because we are caching them.
#
# But set '--no-cache=false' when a picard-private.jar different from
# the default is specified on the assumption that a base image with
# the default jar exists.
#
function runDocker () {
    local -r  jar="$1" timestamp=$(date +"%s")
    local -r tag=$DOCKER_IMAGE_VERSION-$timestamp
    local project=broad-arrays-prod cache=--no-cache=true
    test "$jar" && cache=--no-cache=false
    local -r gcr=us.gcr.io/$project/arrays-picard-private
    echo -e "$gcr:$tag\t$PICARD_PRIVATE_VERSION" >> ../build_arrays_picard_private_docker_version.tsv
    docker build $cache -t $gcr:$tag .
    docker push $gcr:$tag
}

# Run docker login if cannot pull broadinstitute/dsde-toolbox.
#
function ensureLoggedInToDocker () {
    local -r scriptName="$1"
    if docker pull broadinstitute/dsde-toolbox >/dev/null
    then
        : OK
    else
        echo "$scriptName: Are you logged in to Docker?" 1>&2
        docker login
    fi
}

function main () {
    local -r scriptName="${0##*/}" jar="$1"
    local -r tmpdir=$(mktemp -d dockerXXXXXX)
    echo $scriptName: Running: "$0" "$@" 1>&2
    help "$scriptName" "$@"
    echo "Retrieving picard private jar"
    copyPicardPrivateToLocalDockerDirectory "$jar" "$tmpdir"
    cd "$tmpdir"
    makeDockerfile
    ensureLoggedInToDocker "$scriptName"
    runDocker "$jar"
    cd - >/dev/null
    rm -rf "$tmpdir"
}

main "$@"
