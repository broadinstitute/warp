FROM alpine:3.8

ENV TERM=xterm-256color \
    ZCALL_URL=https://github.com/jigold/zCall/raw/master/Version1_AutoCall

ARG ZCALL_VERSION=zCall_Version1.3_AutoCall.zip

WORKDIR /usr/gitc

# Install dependencies.
RUN apk add --no-cache bash python2 unzip curl

# Download zip file and extract contents to working directory
RUN mkdir temp && mkdir zcall && cd temp \
    && curl -LO $ZCALL_URL/$ZCALL_VERSION \
    && unzip $ZCALL_VERSION \
    && cp -R GTC/* /usr/gitc/zcall \
    && rm -rf /usr/gitc/temp
