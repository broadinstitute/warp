FROM us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438

ENV TERM=xterm-256color \
    BOOST_INCLUDEDIR=/usr/include/boost \
    BOOST_LIBRARYDIR=/usr/lib/ \
    HAS_GTEST=0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            bash \
            g++ \
            boost-dev \
    ; \
# Install dragmap
    wget https://github.com/Illumina/DRAGMAP/archive/refs/tags/1.2.1.tar.gz; \
    tar -xf 1.2.1.tar.gz; \
    rm 1.2.1.tar.gz; \
    cd DRAGMAP-1.2.1; \
    make; \
    make install;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
