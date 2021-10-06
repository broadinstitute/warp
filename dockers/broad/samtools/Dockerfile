FROM alpine:3.8

ARG SAMTOOLS_VERSION=1.11

ENV TERM=xterm-256color \
    SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        SAMTOOLS_VERSION=${SAMTOOLS_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            autoconf \
            automake \ 
            bash \
            bzip2-dev \
            curl-dev \
            gcc \
            libressl-dev \
            make \ 
            musl-dev \
            ncurses-dev \
            perl \
            tini \ 
            wget \ 
            xz-dev \
            zlib-dev \
    ; \
# Install samtools
    mkdir temp; \
    cd temp; \
    \
    wget ${SAMTOOLS_URL}; \
    tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2; \
    cd samtools-${SAMTOOLS_VERSION}; \
    \
    ./configure; \
    make; \
    make install; \
    \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
