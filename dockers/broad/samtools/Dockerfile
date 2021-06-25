# Samtools requires lots of Debian only native packages
# This also acts as a base image for many other images that require the JDK
FROM us.gcr.io/broad-dsp-gcr-public/base/jre:8-debian

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

ARG SAMTOOLS_VERSION=1.11

ENV TERM=xterm-256color \
    TINI_VERSION=v0.19.0 \
    SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2

WORKDIR /usr/gitc

# Install dependencies
RUN apt-get update && \
    apt-get -y install wget r-base && \
    \
    # Install samtools
    wget $SAMTOOLS_URL && \
    tar xf samtools-$SAMTOOLS_VERSION.tar.bz2 && \
    rm samtools-$SAMTOOLS_VERSION.tar.bz2 && \
    cd samtools-$SAMTOOLS_VERSION && \
    ./configure && \
    make && \
    make install && \
    cd ../ && \
    rm -r samtools-$SAMTOOLS_VERSION && \
    \
    # Install tini
    wget https://github.com/krallin/tini/releases/download/$TINI_VERSION/tini -O /sbin/tini && \
    chmod +x /sbin/tini && \
    # Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]