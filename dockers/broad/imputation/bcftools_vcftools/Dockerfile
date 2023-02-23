FROM --platform=linux/amd64 python:3.8-alpine

ARG BCFTOOLS_VERSION=1.10.2 \
        VCFTOOLS_VERSION=0.1.16

ENV TERM=xterm-256color

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        BCFTOOLS_VERSION=${BCFTOOLS_VERSION} \
        VCFTOOLS_VERSION=${VCFTOOLS_VERSION}

WORKDIR /usr/gitc

# Install dependencies 
RUN set -eux; \
        apk upgrade; \
        apk add --no-cache \
            autoconf \
            automake \
            bash \
            bzip2-dev \
            curl-dev \
            g++ \
            gcc \
            gsl-dev \
            make \
            musl-dev \
            perl \
            perl-dev \
            tini \
            wget \
            xz-dev \
            zlib-dev \
    ; \
# Install python dependencies
    pip3 install docopt \
    ; \
# Install BCFTools
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2; \
    tar xf bcftools-${BCFTOOLS_VERSION}.tar.bz2; \
    cd bcftools-${BCFTOOLS_VERSION}; \
    \
    ./configure; \
    make; \
    make install; \
    \
    cd ../; \
    rm -r bcftools-${BCFTOOLS_VERSION}; \
    rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    ; \
# Install VCFTools
    wget https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VERSION}/vcftools-${VCFTOOLS_VERSION}.tar.gz; \
    tar xf vcftools-${VCFTOOLS_VERSION}.tar.gz; \
    cd vcftools-${VCFTOOLS_VERSION}; \
    \
    ./configure; \
    make; \
    make install; \
    \
    cd ../; \
    rm -r vcftools-${VCFTOOLS_VERSION}; \
    rm vcftools-${VCFTOOLS_VERSION}.tar.gz \
    ; \
# Install HTSlib
    wget  https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2; \
    tar xf samtools-1.10.tar.bz2; \
    cd samtools-1.10/htslib-1.10; \
    \
    ./configure; \
    make; \
    make install; \
    \
    cd ../..; \
    rm -r samtools-1.10; \
    rm samtools-1.10.tar.bz2

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]


