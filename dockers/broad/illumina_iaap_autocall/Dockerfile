# Using JRE approved image, .NET doesn't exist
FROM us.gcr.io/broad-dsp-gcr-public/base/jre:15-debian

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

ENV TERM=xterm-256color \
    IAAP_CLI_URL=https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/illumina-iaap-cli

ARG IAAP_CLI_VERSION=iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7

WORKDIR /usr/gitc

# Install dependencies
RUN apt-get update && apt-get install -y mono-complete curl && \
    \ 
    # Install Illumina Iaap Cli
    mkdir temp iaap && \
    cd temp && \
    curl -L $IAAP_CLI_URL/$IAAP_CLI_VERSION.tar.gz > iaapcli.tar.gz && \
    tar -xvf iaapcli.tar.gz && \
    cp -R $IAAP_CLI_VERSION/* /usr/gitc/iaap && \
    rm -rf /usr/gitc/temp && \
    \
    # Clear cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*
