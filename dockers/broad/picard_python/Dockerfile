FROM --platform=linux/amd64 python:3.8-alpine

ARG PICARD_PUBLIC_VERSION=2.26.10

ENV TERM=xterm-256color \
    PICARD_URL=https://github.com/broadinstitute/picard/releases/download/${PICARD_PUBLIC_VERSION}/picard.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        PICARD_PUBLIC_VERSION=${PICARD_PUBLIC_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk upgrade; \
        apk add --no-cache \
            bash \
            gcompat \
            libc6-compat \
            openjdk8 \
            tini \
            wget \
    ; \
# Install picard
    wget ${PICARD_URL};

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
