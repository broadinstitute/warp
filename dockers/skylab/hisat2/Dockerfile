FROM --platform=linux/amd64 ubuntu:16.04

ENV TERM=xterm-256color \
        PATH="$PATH:/usr/gitc/hisat2-2.1.0" \
        TINI_VERSION=v0.19.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

RUN set -eux; \
        apt update; \
        apt install -y \
            automake \
            build-essential \
            cmake \
            curl \
            git \
            gcc-multilib \
            liblzma-dev \
            libbz2-dev \
            libboost-all-dev \
            libcurl4-openssl-dev \
            libncurses5-dev \
            libxml2-dev \
            r-base \
            r-base-core \
            r-base-dev \
            unzip \
            wget \
            zip \
            zlib1g-dev \
    ; \
    wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip; \
    unzip hisat2-2.1.0-source.zip; \
    cd hisat2-2.1.0; \
    make; \
    cp hisat2* /usr/local/bin \
    ; \
# Install samtools
    mkdir /usr/local/samtools; \
    cd /usr/local/samtools; \
    wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2; \
    tar -xvf samtools-1.6.tar.bz2; \
    rm samtools-1.6.tar.bz2; \
    cd samtools-1.6; \
    ./configure --prefix=/usr; \
    make; \
    make install \
    ; \
# Install gffread (gtf/gtt tools) \
    mkdir /usr/gitc/gffread; \
    cd /usr/gitc/gffread; \
    git clone https://github.com/gpertea/gclib; \
    git clone https://github.com/gpertea/gffread; \
    cd gffread; \
    make; \
    cp gffread /usr/local/bin/gffread \
    ; \
# Install TINI
    mkdir /usr/gitc/temp; \
    cd /usr/gitc/temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]