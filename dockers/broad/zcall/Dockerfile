# Zcall needs python2
FROM python:2.7.18-alpine3.11

ARG ZCALL_VERSION=zCall_Version1.3_AutoCall

ENV TERM=xterm-256color \
    ZCALL_URL=https://github.com/jigold/zCall/raw/master/Version1_AutoCall

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        ZCALL_VERSION=${ZCALL_VERSION}


WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk --no-cache upgrade; \
        apk add --no-cache \
            bash \
            curl  \
            tini \
            unzip \
    ; \
# Install Zcall
    mkdir temp zcall; \
    cd temp; \
    \
    curl -L $ZCALL_URL/$ZCALL_VERSION.zip > zcall.zip; \
    unzip zcall.zip; \
    cp -R GTC/* /usr/gitc/zcall; \
    \
    rm -rf /usr/gitc/temp

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
