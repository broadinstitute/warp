# Using rocky as temp replacement for CentOS
FROM rockylinux:8.5 

ARG DRAGMAP_VERSION=1.2.1 \
        PICARD_VERSION=2.26.10 \
        SAMTOOLS_VERSION=1.11

ENV TERM=xterm-256color \
        HAS_GTEST=0 \
        TINI_VERSION=v0.19.0 \
        BOOST_LIBRARYDIR=/usr/lib64/boost169 \
        BOOST_INCLUDEDIR=/usr/include/boost169 \
        DRAGMAP_URL=https://github.com/Illumina/DRAGMAP/archive/refs/tags/${DRAGMAP_VERSION}.tar.gz \
        PICARD_URL=https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
        SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        yum upgrade -y; \
        yum install -y epel-release; \
        yum install -y boost169-devel \
            bzip2 \
            bzip2-devel \
            curl-devel \
            gcc \
            gcc-c++ \
            java-1.8.0-openjdk \
            make \
            ncurses-devel \
            openssl-devel \
            wget \
            xz-devel \
            zlib-devel \
        ; \
# Install TINI
    mkdir temp; \
    cd temp; \
    \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini \
    ; \
# Install DRAGMAP
    wget ${DRAGMAP_URL}; \
    tar -xf ${DRAGMAP_VERSION}.tar.gz; \
    rm ${DRAGMAP_VERSION}.tar.gz; \
    \
    cd DRAGMAP-${DRAGMAP_VERSION}; \
    make; \
    mv build/release/dragen-os /usr/local/bin/; \
    \
    cd ../; \
    rm -r DRAGMAP-${DRAGMAP_VERSION} \
    ; \
#Install SAMTOOLS
    wget ${SAMTOOLS_URL}; \
    tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2; \
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2; \
    \
    cd samtools-${SAMTOOLS_VERSION}; \
    make; \
    make install; \
    \
    cd ../; \
    rm -r samtools-${SAMTOOLS_VERSION} \
    ; \
# Install PICARD
    wget ${PICARD_URL}; \
    mkdir /picard; \
    mv picard.jar /picard/;  \
    \ 
    cd ../; \
    rm -r temp \
    ; \
# Clean up cached files
    yum clean all; \
    rm -rf /var/cache/yum;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]