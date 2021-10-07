FROM python:3.8-alpine

ARG BCFTOOLS_VERSION=1.10.2 \
        VCFTOOLS_VERSION=0.1.16

ENV TERM=xterm-256color

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        BCFTOOLS_VERSION=${BCFTOOLS_VERSION} \
        VCFTOOLS_VERSION=${VCFTOOLS_VERSION}

WORKDIR /usr/gitc

# Install dependencies 
RUN set -eux; \
        apk add --no-cache \
            autoconf \
            automake \
            bash \
            bzip2-dev \
            curl-dev \
            g++ \
            gcc \
            gsl-dev \
            libressl-dev \
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
    wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2; \
    tar xf bcftools-1.10.2.tar.bz2; \
    cd bcftools-1.10.2; \
    \
    ./configure; \
    make; \
    make install; \
    \
    cd ../; \
    rm -r bcftools-1.10.2.tar.bz2 \
    ; \
# Install VCFTools
    wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz; \
    tar xf vcftools-0.1.16.tar.gz; \
    cd vcftools-0.1.16; \
    \
    ./configure; \
    make; \
    make install; \
    \
    cd ../; \
    rm -r vcftools-0.1.16.tar.gz 

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]


