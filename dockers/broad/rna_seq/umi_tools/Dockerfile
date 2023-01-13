FROM us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian

ARG UMI_TOOLS_VERSION=1.1.1

ENV TERM=xterm-256color \
        TINI_VERSION=v0.19.0 

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        UMI_TOOLS_VERSION=${UMI_TOOLS_VERSION}

# Install dependencies
RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            gcc \
            wget \
    ; \
# Install umi_tools
    pip3 install umi_tools==${UMI_TOOLS_VERSION} \
    ; \
# Install tini
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini \
    ; \
# Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]