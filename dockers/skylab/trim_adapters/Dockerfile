FROM --platform=linux/amd64 alpine:latest

ARG EAUTILS_VERSION=1.04.807

ENV TERM=xterm-256color \
        EAUTILS_URL=https://github.com/ExpressionAnalysis/ea-utils/archive/refs/tags/${EAUTILS_VERSION}.tar.gz \
        PATH="$PATH:/usr/gitc/ea-utils-${EAUTILS_VERSION}/clipper/"

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        EAUTILS_VERSION=${EAUTILS_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            automake \
            bash \
            build-base \
            cmake \
            gsl-dev \
            tini \
            wget \
            zlib-dev \
        ; \
# Install ea-utils  \
    wget ${EAUTILS_URL}; \
    tar -xzvf ${EAUTILS_VERSION}.tar.gz; \
    cd ea-utils-${EAUTILS_VERSION}/clipper; \
    make; \
    make install;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--"]