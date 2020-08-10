#!/bin/bash
set -e

BASE_DIR=$(cd $(dirname $0) && pwd)
TMPDIR=$(mktemp -d dockerXXXXXX)

# Update DOCKER_IMAGE_VERSION after any substantial changes to the
# image this builds.
#
declare -r DOCKER_IMAGE_VERSION=1.0.1

# Update this when there is a new release of Picard to use as the
# default jar.
#
declare -r PICARD_PRIVATE_VERSION=1.1413

# I don't know why vpicard05:/..., but these are scp sources below.
#
declare -r DEFAULT_SOFTWARE=vpicard05:/seq/software
declare -r DEFAULT_PICARD_PRIVATE_JAR="$DEFAULT_SOFTWARE/picard/${PICARD_PRIVATE_VERSION}/bin/picard-private.jar"

# Show usage and fail if there is a problem with the command line.
#
function help {
    local scriptName="$1" environment="$2" jar="$3" ok=no
    local -r -a environments=(dev staging prod)
    local -r -a tools=(docker gcloud vault wget)
    local -r -a usage=(
        ''
        "$scriptName: Build a GotC Multi Sample Arrays docker image and push it to GCR."
        ''
        "Usage: $scriptName <environment> [<jar>]"
        ''
        "Where: <environment> is one of: ${environments[*]}"
        '       <jar> is the path to a picard-private.jar.'
        ''
        "The default <jar> is $DEFAULT_PICARD_PRIVATE_JAR"
        ''
        "Example: $scriptName dev ./dist/picard-private.jar"
        ''
        "Note: You need at least these tools installed: ${tools[*]}"
        ''
    )
    for e in ${environments[*]}; do test $e = "$environment" && ok=yes; done
    test $ok = no && echo "${scriptName}: Unrecognized environment: '$environment'" 1>&2
    for t in "${tools[@]}"; do which $t >/dev/null || ok=no; done
    if test $ok = no
    then
        for line in "${usage[@]}"; do echo "$line" 1>&2; done
        exit 1
    fi
}

# Copy the (possibly remote) picard-private.jar binaries
# to the docker build directory.
#
function copyPicardPrivateToDockerBuildDirectory () {
    local -r jar="$1"
    local picard="$jar"
    test "$picard" || picard="$DEFAULT_PICARD_PRIVATE_JAR"
    scp "$picard" ${TMPDIR}/
}

# Pull the credentials for the gcr-arrays-reader-only from the vault.
# Note that we do the docker business here in order to run 'jq' which
# is used to format the json properly.
#
function runDockerToPullCredentialsFromVault () {
    local -r environment=$1
    local -r working="${PWD}:/working" root="${HOME}:/root"
    local -r dev=broadinstitute/dsde-toolbox:dev
    local -r account=gcr-arrays-reader-only-account.json
    local -r source=secret/dsde/gotc/$environment/gcr/$account
    local -r destination=./${TMPDIR}/$account
    local -r cmd="vault read --format=json $source | jq .data > $destination"
    docker run -it --rm -v "$working" -v "$root" $dev /bin/bash -c "$cmd"
}

# Make the Dockerfile.
#
function makeDockerfile () {
    local -r environment=$1
    echo "FROM centos:centos6

ENV TERM=xterm-256color

WORKDIR /usr/gitc

RUN rpm -Uvh http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm \
  && yum update -y \
  && yum install -y curl which wget tar glib2-devel cairo-devel libtiff-devel libjpeg-devel giflib-devel libexif-devel libX11 libX11-devel libXt libXt-devel java-1.8.0-openjdk R \
  && yum groupinstall -y 'Development Tools' \
  && yum install -y zlib-devel openssl-devel sqlite-devel bzip2-devel

# put python 2.7.10 on this image so gcloud will work (deprecated as of January 2017)
RUN curl -O https://www.python.org/ftp/python/2.7.10/Python-2.7.10.tgz \
    && tar xvfz Python-2.7.10.tgz \
    && cd Python-2.7.10 \
    && ./configure --enable-shared --prefix=/usr/local LDFLAGS=-Wl,--rpath=/usr/local/lib \
    && make \
    && make altinstall \
    && cd .. \
    && rm -rf Python-2.7.10 \
    && rm Python-2.7.10.tgz

# tell cloudsdk which python to use since we dont want to mess with the system python
ENV CLOUDSDK_PYTHON=/usr/local/bin/python2.7

# download the gsutil install script
RUN curl https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash > gsutil_install.sh

# install gsutil
RUN /bin/bash gsutil_install.sh --disable-prompts


# restart shell
RUN exec -l /bin/bash

COPY . .

# activate public service account
RUN /root/google-cloud-sdk/bin/gcloud auth activate-service-account gcr-arrays-reader-only@broad-gotc-$environment.iam.gserviceaccount.com --key-file gcr-arrays-reader-only-account.json
" > ./Dockerfile
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
    local -r environment=$1 jar="$2" timestamp=$(date +"%s")
    local -r tag=$environment-$DOCKER_IMAGE_VERSION-$timestamp
    local project=broad-arrays cache=--no-cache=true
    test $environment = dev && project=broad-gotc
    project=$project-$environment
    test "$jar" && cache=--no-cache=false
    local -r gcr=us.gcr.io/$project/multi-sample-arrays
    echo -e "$tag\t$PICARD_PRIVATE_VERSION" >> ../build_docker_version.tsv
    docker build $cache -t broadinstitute/multi-sample-arrays:$tag .
    docker tag broadinstitute/multi-sample-arrays:$tag $gcr:$tag
    gcloud docker -- push $gcr:$tag
    # Replace docker image in options file
    cd ..
    sed -i '' "s/\"docker\":\ \"[^\"]*\"/\"docker\":\ \"us.gcr.io\/$project\/multi-sample-arrays:$tag\"/g" MultiSampleArrays.$environment.options.json
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
    local -r scriptName="${0##*/}" environment="$1" jar="$2"
    echo $scriptName: Running: "$0" "$@" 1>&2
    help "$scriptName" "$@"
    ensureLoggedInToDocker "$scriptName"
    copyPicardPrivateToDockerBuildDirectory "$jar"
    runDockerToPullCredentialsFromVault "$environment"
    cd ./${TMPDIR}
    makeDockerfile "$environment"
    runDocker "$environment" "$jar"
    cd ${BASE_DIR}
    rm -rf ${TMPDIR}
}

main "$@"
