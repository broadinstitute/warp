FROM --platform=linux/amd64 us.gcr.io/broad-gotc-prod/snaptools:1.0.0-1.4.8-1660834594

ARG BWA_VERSION=0.7.17

ENV TERM=xterm-256color \
       BWA_URL=https://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2/download

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

RUN set -eux; \
    apt-get update; \
    apt-get upgrade -y; \
    mkdir temp; \
    cd temp; \
    wget -O "bwa.tar.bz2" ${BWA_URL}; \
    tar -xvjf bwa.tar.bz2; \
    cd bwa-${BWA_VERSION}; \
    make; \
    cp bwa /usr/local/bin;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]