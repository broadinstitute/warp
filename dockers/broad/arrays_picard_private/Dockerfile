FROM frolvlad/alpine-glibc:glibc-2.33

ARG PICARD_PRIVATE_VERSION=257537c72dae29257b09bacf413505eed295ac32

ENV TERM=xterm-256color \
    NO_VAULT=true \
    ARTIFACTORY_URL=https://broadinstitute.jfrog.io/artifactory/libs-snapshot-local/org/broadinstitute/picard-private

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        PICARD_PRIVATE_VERSION=${PICARD_PRIVATE_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk --no-cache upgrade; \
        apk add --no-cache \    
            bash \
            curl \ 
            findutils \
            jq \ 
            openjdk8-jre \
            python3 \
            tini \
            unzip \
            wget \
    ; \
# Download picard private jar
    curl -L ${ARTIFACTORY_URL}/${PICARD_PRIVATE_VERSION}/jars/picard-private-all-${PICARD_PRIVATE_VERSION}.jar > picard-private.jar \
    ; \
# Download the gsutil install script
    wget https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash -O - | bash \
    ; \ 
# Set up Vault
    curl -L https://releases.hashicorp.com/vault/1.0.2/vault_1.0.2_linux_amd64.zip > temp.zip; \
    unzip temp.zip; \
    rm temp.zip; \ 
    mv vault /usr/local/bin/

# Set tini as default entry point
ENTRYPOINT ["/sbin/tini", "--"]
