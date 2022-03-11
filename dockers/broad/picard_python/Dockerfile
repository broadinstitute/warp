FROM python:3.8-alpine

ARG PICARD_PUBLIC_VERSION=2.26.10

ENV TERM=xterm-256color \
    PICARD_URL=https://github.com/broadinstitute/picard/releases/download/${PICARD_PUBLIC_VERSION}/picard.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        PICARD_PUBLIC_VERSION=${PICARD_PUBLIC_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            gcompat \
            libc6-compat \
            openjdk8 \
    ; \
# Install picard
    wget ${PICARD_URL};

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
