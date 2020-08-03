#!/bin/bash
set -e

# Update DOCKER_IMAGE_VERSION after any substantial changes to the
# image this builds.
#
declare -r DOCKER_IMAGE_VERSION=4.0.1

# Show usage and fail if there is a problem with the command line.
#
function help {
    local scriptName="$1" ok=yes
    local -r -a tools=(docker wget jq)
    local -r -a usage=(
        ''
        "$scriptName: Build a Zcall docker image and push it to the GCR."
        ''
        'The GCR is the Google Cloud Registry.'
        ''
        "Usage: $scriptName"
        ''
        "Example: $scriptName"
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

# Get zCall from Github into the local docker directory.
#
function getZcallFromGitHub () {
    local -r zcallzip=zCall_Version1.3_AutoCall.zip
    wget https://github.com/jigold/zCall/raw/master/Version1_AutoCall/$zcallzip
    unzip $zcallzip
    rm $zcallzip
    mv ./GTC/* ./zcall
    rmdir GTC
}

# Make the Dockerfile.
#
function makeDockerfile () {
cat << EOF > ./Dockerfile
FROM alpine:3.8

ENV TERM=xterm-256color

WORKDIR /usr/gitc

# Install dependencies.
RUN apk --update add bash python2

COPY . .
EOF
}

# Run docker build, then tag and push the image to the GCR.
#
# By default, set `--no-cache=true` here so we will always build our
# images from scratch instead of possibly missing failing steps in the
# process because we are caching them.
#
function runDocker () {
    local -r timestamp=$(date +"%s")
    local -r tag=$DOCKER_IMAGE_VERSION-$timestamp
    local project=broad-gotc-prod cache=--no-cache=true
    local -r gcr=us.gcr.io/$project/zcall
    echo -e "$gcr:$tag" >> ../build_zcall_docker_version.tsv
    docker build $cache -t $gcr:$tag .
    gcloud docker -- push $gcr:$tag
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
    local -r scriptName="${0##*/}"
    local -r tmpdir=$(mktemp -d dockerXXXXXX)
    echo $scriptName: Running: "$0" "$@" 1>&2
    help "$scriptName" "$@"
    mkdir "$tmpdir/zcall"
    cd "$tmpdir"
    getZcallFromGitHub
    makeDockerfile
    ensureLoggedInToDocker "$scriptName"
    runDocker
    cd - >/dev/null
    rm -rf "$tmpdir"
}

main "$@"
