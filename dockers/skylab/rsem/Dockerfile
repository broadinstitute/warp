FROM --platform=linux/amd64 ubuntu:16.04

ENV TERM=xterm-256color \
        RSEM_VERSION=1.3.0 \
        STAR_VERSION=2.5.3a \
        TINI_VERSION=v0.19.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

# Install compiler, perl , R and stuff
RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            apt-utils \
            automake \
            bowtie \
            bowtie2 \
            build-essential \
            cmake \
            curl \
            gcc-multilib \
            libboost-all-dev \
            libbz2-dev \
            libncurses5-dev \
            libxml2-dev \
            perl \
            perl-base \
            r-base \
            r-base-core \
            r-base-dev \
            wget \
            zlib1g-dev \
    ; \
# Install STAR \
    mkdir /usr/gitc/temp; \
    cd /usr/gitc/temp; \
    wget https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz; \
    tar -xf ${STAR_VERSION}.tar.gz; \
    cp STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin \
    ; \
# Install RSEM \
    cd /usr/gitc; \
    wget https://github.com/deweylab/RSEM/archive/refs/tags/v${RSEM_VERSION}.tar.gz; \
    tar -xf v${RSEM_VERSION}.tar.gz; \
    cd RSEM-${RSEM_VERSION}; \
    make; \
    make ebseq; \
    make install \
    ; \
# Install TINI
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]