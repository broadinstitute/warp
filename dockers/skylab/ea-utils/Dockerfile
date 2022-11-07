FROM --platform=linux/amd64 ubuntu:16.04

ARG EAUTILS_VERSION=1.04.807

ENV TERM=xterm-256color \
        EAUTILS_URL=https://github.com/ExpressionAnalysis/ea-utils/archive/refs/tags/${EAUTILS_VERSION}.tar.gz \
        EAUTILS_VERSION=${EAUTILS_VERSION} \
        TINI_VERSION=v0.19.0 \
        PATH="$PATH:/usr/gitc/ea-utils-${EAUTILS_VERSION}/clipper/"

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        EAUTILS_VERSION=${EAUTILS_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
    apt update; \
    apt install -y \
       automake \
       build-essential \
       cmake \
       libgsl0-dev \
       wget \
       zlib1g-dev \
    ; \
# Install ea-utils
    wget ${EAUTILS_URL}; \
    tar -xzvf ${EAUTILS_VERSION}.tar.gz; \
    cd ea-utils-${EAUTILS_VERSION}/clipper; \
    make; \
    make install \
    ; \
# Install tini
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]