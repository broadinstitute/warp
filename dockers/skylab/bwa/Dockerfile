FROM --platform=linux/amd64 ubuntu:16.04

ARG BWA_VERSION=0.7.17

ENV TERM=xterm-256color \
        TINI_VERSION=v0.19.0 \
        BWA_URL=https://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2/download

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        BWA_VERSION=${BWA_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
        apt update; \
        apt install -y \
            build-essential \
            bzip2 \
            zlib1g \
            zlib1g-dev \
            wget \
    ; \
# Install BWA \
    mkdir temp; \
    cd temp; \
    wget -O "bwa.tar.bz2" ${BWA_URL}; \
    tar -xvjf bwa.tar.bz2; \
    cd bwa-${BWA_VERSION}; \
    make; \
    cp bwa /usr/local/bin \
    ; \
# Install TINI
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]
