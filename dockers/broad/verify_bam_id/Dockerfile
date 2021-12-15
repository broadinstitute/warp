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
            libcurl4-gnutls-dev \
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
