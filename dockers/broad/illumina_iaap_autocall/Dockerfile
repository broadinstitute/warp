FROM frolvlad/alpine-mono:5.4-glibc

ARG IAAP_CLI_VERSION=iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7

ENV TERM=xterm-256color \
    IAAP_CLI_URL=https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/illumina-iaap-cli

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        IAAP_CLI_VERSION=${IAAP_CLI_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk --no-cache upgrade; \
        apk add --no-cache \ 
            curl \
            bash \
            icu \
            libstdc++ \
            tini \
    ; \ 
# Install Illumina Iaap Cli
    mkdir temp iaap; \
    cd temp; \
    \
    curl -L ${IAAP_CLI_URL}/${IAAP_CLI_VERSION}.tar.gz > iaapcli.tar.gz; \
    tar -xf iaapcli.tar.gz; \
    cp -R ${IAAP_CLI_VERSION}/* /usr/gitc/iaap; \
    \
    rm -rf /usr/gitc/temp

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]
