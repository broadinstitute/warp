FROM frolvlad/alpine-glibc

ENV TERM=xterm-256color

WORKDIR /usr/gitc

# Install dependencies.
RUN apk --update add bash curl findutils jq wget openjdk8-jre python3 unzip

# download the gsutil install script
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash -O - | bash

# Set up Vault
ENV NO_VAULT=true
RUN curl https://releases.hashicorp.com/vault/1.0.2/vault_1.0.2_linux_amd64.zip > temp.zip \
    && unzip temp.zip \
    && rm temp.zip \
    && mv vault /usr/local/bin/

COPY . .
