# Eagle needs to run on debian based system
FROM us.gcr.io/broad-dsp-gcr-public/base/python:debian

ARG EAGLE_VERSION=2.4

ENV TERM=xterm-256color \
    TINI_VERSION=v0.19.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        EAGLE_VERSION=${EAGLE_VERSION} 
        

WORKDIR /usr/gitc

# Install dependencies 
RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            build-essential \
            libgomp1 \
            libbz2-dev \
            liblzma-dev  \
            libz-dev \
            pkg-config \
            tabix \
            wget \
    ; \
# Install Eagle
    wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/old/Eagle_v${EAGLE_VERSION}.tar.gz; \
    tar xf Eagle_v${EAGLE_VERSION}.tar.gz; \
    \
    mv Eagle_v${EAGLE_VERSION}/eagle .; \
    rm -r Eagle_v${EAGLE_VERSION}; \
    rm Eagle_v${EAGLE_VERSION}.tar.gz; \
    \
    chmod +x eagle \
    ; \
# Download genetic map
    wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz \
    ; \
# Install tini
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini \
    ; \
# Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
