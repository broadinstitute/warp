FROM alpine:3.8

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

ENV TERM=xterm-256color \
    ZCALL_URL=https://github.com/jigold/zCall/raw/master/Version1_AutoCall

ARG ZCALL_VERSION=zCall_Version1.3_AutoCall

WORKDIR /usr/gitc

# Install dependencies
RUN apk add --no-cache bash python2 unzip curl && \
    \
    # Download zip file and extract contents to working directory
    mkdir temp zcall && cd temp && \
    curl -L $ZCALL_URL/$ZCALL_VERSION.zip > zcall.zip  && \
    unzip zcall.zip  && \
    cp -R GTC/* /usr/gitc/zcall && \
    rm -rf /usr/gitc/temp
