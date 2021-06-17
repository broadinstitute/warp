FROM frolvlad/alpine-glibc

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

ENV TERM=xterm-256color \
    NO_VAULT=true

WORKDIR /usr/gitc

# Install dependencies.
RUN apk add --no-cache bash curl findutils jq wget openjdk8-jre python2 unzip && \
    \
    # Download the gsutil install script
    wget https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash -O - | bash && \
    \ 
    # Set up Vault
    curl -L https://releases.hashicorp.com/vault/1.0.2/vault_1.0.2_linux_amd64.zip > temp.zip && \
    unzip temp.zip && \
    rm temp.zip  && \ 
    mv vault /usr/local/bin/

COPY picard-private.jar .
