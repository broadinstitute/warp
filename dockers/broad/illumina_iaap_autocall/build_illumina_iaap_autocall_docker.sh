#!/bin/bash
set -e

declare -r DOCKER_IMAGE_VERSION=1.0.2

# rsync iaap autocall dlls
#
declare -r DEFAULT_SOFTWARE=picard03:/seq/software/iaap/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7
declare -r DEFAULT_IAAP_BINS="$DEFAULT_SOFTWARE/*"

# Show usage and fail if there is a problem with the command line.
#
function help {
    local scriptName="$1" ok=yes
    local -r -a tools=(docker jq)
    local -r -a usage=(
        ''
        "$scriptName: Build an Illumina IAAP Autocall docker image and push it to the GCR."
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

# Copy the AutoCall binaries to the local docker build directory ($dir).
# Don't copy AutoCall binaries if there are already a bunch of files in (docker
# build)/iaap/ ... which is a terrible hack, BTW.
#
function copyIaapAutocallToLocalDockerDirectory () {
    local -r dir="$1"
    local -r -a files=("$dir"/iaap/*)
    let "${#files[*]} > 40" || rsync -a --progress "$DEFAULT_IAAP_BINS" "$dir"/iaap
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
    local -r  timestamp=$(date -u +"%Y-%m-%d_%H-%M-%SZ")
    local -r tag=$DOCKER_IMAGE_VERSION-$timestamp
    local project=broad-gotc-prod cache=--no-cache=true
    local -r gcr=us.gcr.io/$project/illumina-iaap-autocall
    docker build $cache -t $gcr:$tag .
    docker push $gcr:$tag
    echo -e "$gcr:$tag" >> ../build_illumina_iaap_autocall_docker_version.tsv
}

function main () {
    local -r scriptName="${0##*/}"
    local -r tmpdir=$(mktemp -d dockerXXXXXX)
    echo $scriptName: Running: "$0" "$@" 1>&2
    help "$scriptName" "$@"
    mkdir "$tmpdir/iaap"
    copyIaapAutocallToLocalDockerDirectory "$tmpdir"
    cp Dockerfile "$tmpdir"
    cd "$tmpdir"
    # ensureLoggedInToDocker "$scriptName"
    runDocker
    cd - >/dev/null
    rm -rf "$tmpdir"
}

main "$@"
