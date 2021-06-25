FROM alpine:3.8

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

ARG SAMTOOLS_VERSION=1.11

ENV TERM=xterm-256color \
    SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2

WORKDIR /usr/gitc

# Install dependencies
RUN apk add --no-cache \
    wget autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev libressl-dev ncurses-dev tini && \
    \
    # Install samtools
    mkdir temp && \
    cd temp && \
    wget $SAMTOOLS_URL && \
    tar xf samtools-$SAMTOOLS_VERSION.tar.bz2 && \
    cd samtools-$SAMTOOLS_VERSION && \
    ./configure && \
    make && \
    make install && \
    rm -r /usr/gitc/temp

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]