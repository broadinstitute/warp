FROM python:3.7
MAINTAINER Broad Institute DSDE <dsde-engineering@broadinstitute.org>

WORKDIR /usr/gitc

ARG GIT_HASH

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y zlib1g-dev libcurl4-openssl-dev libssl-dev libbz2-dev libhts-dev gcc g++ cmake autoconf wget unzip git

RUN wget https://github.com/Griffan/VerifyBamID/archive/$GIT_HASH.zip && \
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
    rm -rf $GIT_HASH.zip VerifyBamID-$GIT_HASH
