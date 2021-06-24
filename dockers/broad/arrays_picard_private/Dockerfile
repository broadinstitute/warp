FROM frolvlad/alpine-glibc:glibc-2.33

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

ENV TERM=xterm-256color \
    NO_VAULT=true \
    ARTIFACTORY_URL=https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/picard-private

ARG PICARD_PRIVATE_VERSION=61af9bff4587783e5981a496f422ea36102482b5

WORKDIR /usr/gitc

# Install dependencies
RUN apk add --no-cache bash curl findutils jq wget tini openjdk8-jre python2 unzip && \
    \
    # Download picard private jar
    curl -L $ARTIFACTORY_URL/$PICARD_PRIVATE_VERSION/picard-private-$PICARD_PRIVATE_VERSION.jar > picard-private.jar && \
    \
    # Download the gsutil install script
    wget https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash -O - | bash && \
    \ 
    # Set up Vault
    curl -L https://releases.hashicorp.com/vault/1.0.2/vault_1.0.2_linux_amd64.zip > temp.zip && \
    unzip temp.zip && \
    rm temp.zip  && \ 
    mv vault /usr/local/bin/

# Set tini as default entry point
ENTRYPOINT ["/sbin/tini", "--"]