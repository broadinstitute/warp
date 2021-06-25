FROM us.gcr.io/broad-dsp-gcr-public/base/python:debian

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

ENV TERM=xterm-256color \
    BAMID_URL=https://github.com/Griffan/VerifyBamID/archive

WORKDIR /usr/gitc

ARG GIT_HASH=c1cba76e979904eb69c31520a0d7f5be63c72253

# Install dependencies
RUN apt-get update && \
<<<<<<< HEAD:dockers/broad/VerifyBamId/Dockerfile
    apt-get upgrade -y && \
    apt-get install -y zlib1g-dev libcurl4-openssl-dev libssl-dev libbz2-dev libhts-dev gcc g++ cmake autoconf wget unzip git  && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/Griffan/VerifyBamID/archive/$GIT_HASH.zip && \
=======
    apt-get install -y zlib1g-dev libcurl4-openssl-dev libssl-dev libbz2-dev libhts-dev gcc g++ cmake autoconf wget unzip git && \
    \
    # Install BamID
    wget $BAMID_URL/$GIT_HASH.zip && \
>>>>>>> trivy-action:dockers/broad/verify_bam_id/Dockerfile
    unzip $GIT_HASH.zip && \
    cd VerifyBamID-$GIT_HASH && \
    mkdir build && \
    cd build && \
    CC=$(which gcc) CXX=$(which g++) cmake ..  && \
    cmake && \
    make && \
    make test && \
    cd ../../ && \
    mv VerifyBamID-$GIT_HASH/bin/VerifyBamID . && \
    rm -rf $GIT_HASH.zip VerifyBamID-$GIT_HASH && \
    \
    # Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*
