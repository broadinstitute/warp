FROM debian:bullseye-slim

ARG FASTP_VERSION=0.20.1

ENV TERM=xterm-256color \
        FASTP_URL=http://opengene.org/fastp/fastp.${FASTP_VERSION} \
        TINI_VERSION=v0.19.0 

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        FASTP_VERSION=${FASTP_VERSION}

WORKDIR /usr/gitc

# Install dependencies

RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            wget \
    ; \
# Install fastp
    wget ${FASTP_URL} ;\
    \
    mv fastp.${FASTP_VERSION} /usr/local/bin/fastp; \
    chmod a+x /usr/local/bin/fastp \
    ; \
# Install tini 
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini \
    ; \
# Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]

