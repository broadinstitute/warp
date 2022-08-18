FROM --platform=linux/amd64 python:3.7.2

ARG BWA_VERSION=0.7.17

ENV TERM=xterm-256color \
        BWA_URL=https://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2/download \
        TINI_VERSION=v0.19.0 \
        PATH=$PATH:/usr/gitc

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

RUN set -eux; \
    pip install --upgrade pip; \
    apt-get update; \
    apt-get upgrade -y \
    ; \
    # Install snaptools
    pip install snaptools; \
    mkdir -p /home/r3fang/anaconda2/bin/; \
    ln -s `which python` /home/r3fang/anaconda2/bin/python \
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