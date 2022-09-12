# Docker Style Guide 

This style guide provides formatting guidelines and best practices for writing Dockerfiles in WARP.

## :book: Table of Contents

* [Overview](#overview)  
* [Goals](#goals)
  * [Small images](#small) 
    * [Alpine base](#alpine)
    * [Specifying image platform](#platform)
    * [Minimal RUN steps](#minimal-run)
  * [Publicly accessible](#publicly)
  * [Image scanning](#scanning)
  * [Semantic tagging](#semantic)
  * [Proper process reaping](#process)
* [Build Scripts and README](#build)
* [Formatting](#formatting)
* [Troubleshooting and running standalone](#trouble)
## <a name="overview"></a> Overview

WARP maintains a collection of docker images which are used as execution environments for various cloud-optimized data processing pipelines. Many of these image require specific sets of tools and dependencies to run and can be thought of as _custom_ images rather than traditional application images. 
Building and maintaining these images can be challenging; this document provides a set of guidelines to assist in developing docker images for WARP.

## <a name="goals"></a> Goals

The following are some goals/guidelines we want to strive for when writing our Dockerfiles.

### <a name="small"></a> Small images

Building a smaller image offers advantages such as faster upload and download times along with reduced storage costs and minimized attack vector. Two of the easiest ways to minimize the size of your image is to use a small base image and to reduce the number of layers in your image.

#### <a name="alpine"></a> Alpine base

The easiest way to have a small image is to use an [Alpine](https://alpinelinux.org/) base image. Alpine linux, compared to Debian, RHEL etc., is designed specifically for security and resource efficiency, and lends itself perfectly to be used as a building block for Docker images.

Along with being a small base, Alpine also has built in deletion of package index and provides [tini](https://github.com/krallin/tini) natively through APK.

There are some instances where a Debian base image is unavoidable, specifically in the case where dependencies don't exist in APK. It is suggested that you only go to a Debian base as a last resort.


##### :eyes: Example
```dockerfile

#OKAY, NOT GREAT - uses debian
FROM python:debian

RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            curl \
            bash \
    ; \
# Must clean up cache manually with Debian
    apt-get clean && rm -rf /var/lib/apt/list/*

# GOOD - uses alpine
FROM alpine:3.9

RUN set -eux; \
        apk add --no-cache \
            curl \
            bash \
```

#### <a name="platform"></a> Specifying image platform

Docker images built on ARM-based machines such as the new M-series Macs may run into execution issues with our automated PR test suite.
One way to avoid these issues is to use a `linux/amd64` base image by including the `--platform="linux/amd64` flag after the `FROM` keyword.

##### :eyes: Example
```dockerfile
# Use the amd64 version of alpine
FROM --platform="linux/amd64" alpine
```

#### <a name="minimal-run"></a> Minimal RUN steps

Having minimal `RUN`steps (ideally one) is another highly effective way to reduce the size of your image. Each instruction in a Dockerfile creates a [layer](https://docs.docker.com/storage/storagedriver/) and these layers are what add up to build the final image.
When you use multiple `RUN` steps it creates additional unnecessary layers and bloats your image.

An alternative to having a single `RUN` step is to use [multi-stage builds](https://docs.docker.com/develop/develop-images/multistage-build/) which are effective when the application you are containerizing is just a statically linked binary. 
Just to note, many of the images maintained in WARP require a handful of system-level dependencies and custom packages so multi-stages builds are typically not used.

##### :eyes: Example
```dockerfile

# BAD - uses multiple RUN steps
RUN set -eux
RUN apk add --no-cache curl bash wget
RUN wget https://www.somezipfile.com/zip
RUN unzip zip

# GOOD - uses single RUN step
RUN set -eux; \
        apk add --no-cache \
            curl \
            bash \
    ; \
    wget https://www.somezipfile.com/zip; \
    unzip zip
```

### <a name="publicly"></a> Publicly accessible

The pipelines that we maintain in WARP are designed for public use, ideally we would like our docker images to be publicly available as well. This would mean the following conditions must be true.

* Anybody can pull our images
* Anybody can build our images

For anybody to be able to pull our images they must be hosted on a public container registry, we host all of our images in public repos on GCR (our 'official' location) and Quay (for discoverability).

* GCR - `us.gcr.io/broad-gotc-prod`
* Quay - `quay.io/broadinstitute/broad-gotc-prod`

For anybody to be able to build our images, all functionality should be encapsulated in the Dockerfile. Any custom software packages, dependencies etc. have to be downloaded from public links within the Dockerfile, this obviously means that we should not be copying files from within the Broad network infrastructure into our images.

### <a name="scanning"></a> Image scanning


All images that we build are scanned for critical vulnerabilities on every pull request. For this we use a github-action that leverages [trivy](https://github.com/aquasecurity/trivy) for scanning. If you build a new image please add it to the action [here](../.github/workflows/trivy.yml).

### <a name="semantic"></a> Semantic tagging


We recommend against using rolling tags like `master` or `latest` when building images. Rolling tags make it hard to track down versions of images since the underlying image hash and content could be different across the same tags. Instead, we ask that you use a semantic tag that follows the convention below:

##### `us.gcr.io/broad-gotc-prod/samtools:<image-version>-<samtools-version>-<unix-timestamp>` 

This example is for an image we use containing `samtools`. The 'image-version' in this case is the traditional `major.minor.patch` version of the image being built, which is updated when changes to the image (underlying OS, system level packages, etc.) unrelated to `samtools` are made. The 'samtools-version' here correlates with the specific version of `samtools` being used, having this information in the tag makes it easy for users to identify and not have to track down. Lastly, a unix timestamp in included to avoid any potential issues with Cromwell image caching.

### <a name="process"></a> Proper process reaping


Classic init systems like systemd are used to reap orphaned, zombie processes. Typically, these orphaned processes are reattached to the process at PID 1 which will reap them when they die. In a container this responsibility falls to process at PID 1 which is by default `/bin/sh`...this obviously will not handle process reaping. Because of this you run the risk of expending excess memory or resources within your container. A simple solution to this is to use `tini` in all of our images, a lengthy explanation of what this package does can be found [here](https://github.com/krallin/tini/issues/8).

Luckily `tini` is available natively through APK so all you have to do is install it and set it as the default entrypoint!

##### :eyes: Example
```dockerfile

FROM alpine:3.9

RUN set -eux; \
        apk add --no-cache \
            tini

ENTRYPOINT ["/sbin/tini" , "--"]
```

## <a name="build"></a> Build Scripts and README

To make life easier when building and pushing our images we like to have an easy-to-use `docker_build.sh` that sits next to each Dockerfile. These scripts should have configurable inputs for the version of tools (samtools, picard, zcall, etc.) being used in the image. Additionally, we like to keep a record of the versions built and being used by writing the images to the accompanying `docker_versions.tsv`, this should be done automatically by your build script.

For first time users of these images it is helpful to have a README which gives a high-level overview of the image.

See the examples for samtools([docker_build](./broad/samtools/docker_build.sh), [docker_versions](./broad/samtools/docker_versions.tsv), [README](./broad/samtools/README.md))

## Formatting

Formatting our Dockerfiles consistently helps improve readability and eases maintenance headaches down the road. The following are a couple of tenants that we follow when writing our Dockerfiles:

* ARGS, ENV, LABEL in that order
* Always add versions of tools in the LABEL
* Single RUN steps
* Alphabetize package install
* Clean up package index cache
* Use ; instead of && for line continuation
* Logically separate steps within RUN
* Four spaces per tab indent
* Short comments to describe each step
* tini is always default entrypoint

The following is a good example for our `verify_bam_id` image. This Dockerfile shows how to install packages, install tini and clean up cached index files for a debian base image.

##### :eyes: Example
```dockerfile
# Have to use debian based image, many of the installed packages here are not available in Alpine
FROM us.gcr.io/broad-dsp-gcr-public/base/python:debian

ARG GIT_HASH=c1cba76e979904eb69c31520a0d7f5be63c72253

ENV TERM=xterm-256color \
    BAMID_URL=https://github.com/Griffan/VerifyBamID/archive \
    TINI_VERSION=v0.19.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        GIT_HASH=${GIT_HASH}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            autoconf \
            cmake \
            g++ \
            gcc \
            git \
            libbz2-dev \
            libcurl4-openssl-dev \
            libhts-dev \
            libssl-dev  \
            unzip \
            wget \
            zlib1g-dev \
    ; \
# Install BamID
    wget ${BAMID_URL}/${GIT_HASH}.zip; \
    unzip ${GIT_HASH}.zip; \
    \
    cd VerifyBamID-${GIT_HASH}; \
    mkdir build;  \
    cd build; \
    CC=$(which gcc) CXX=$(which g++) cmake ..; \
    \
    cmake; \
    make; \
    make test; \
    \
    cd ../../; \
    mv VerifyBamID-${GIT_HASH}/bin/VerifyBamID .; \
    rm -rf ${GIT_HASH}.zip VerifyBamID-${GIT_HASH} \
    ; \
# Install tini
    wget https://github.com/krallin/tini/releases/download/$TINI_VERSION/tini -O /sbin/tini; \
    chmod +x /sbin/tini \
    ; \
# Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
```

## <a link="trouble"></a> Troubleshooting and running standalone

The WARP dockers are designed to be run from their respective WDL pipelines. However, if you need to run a Docker independent of a WDL for testing or troubleshooting, you'll likely need to explicity instruct it to run a `bash` shell in the `run` command. An example of this is shown in the terminal command below: 

```bash
docker run -it --rm <docker url> bash
```

If you have any questions or would like some more guidance on writing Dockerfiles please file a [GitHub issue in WARP](https://github.com/broadinstitute/warp/issues/new).
