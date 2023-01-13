FROM adoptopenjdk/openjdk8:alpine-slim

# TODO update to 4.2.6.1 and build
ARG GATK4_VERSION=4.2.6.1
ARG GATK3_VERSION=3.5

ENV TERM=xterm-256color \
    GATK4_URL=https://github.com/broadinstitute/gatk/releases/download/${GATK4_VERSION}/gatk-${GATK4_VERSION}.zip \
    GATK3_URL=https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/gatk/gatk-utils/${GATK3_VERSION}/gatk-utils-${GATK3_VERSION}.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        GATK4_VERSION=${GATK4_VERSION} \
        GATK3_VERSION=${GATK3_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            bash \
            python3 \
            py3-pip \
            tini \
            unzip \
            wget \
    ; \
# Install GATK
    mkdir gatk4; \
    mkdir temp; \
    cd temp; \
    \
    # gatk3
    wget ${GATK3_URL}; \
    cp gatk-utils-${GATK3_VERSION}.jar GATK35.jar; \
    mv GATK35.jar /usr/gitc; \
    \
    # gatk4
    wget ${GATK4_URL}; \
    unzip gatk-${GATK4_VERSION}.zip; \
    cp -r gatk-${GATK4_VERSION}/* /usr/gitc/gatk4/; \
    chmod +x /usr/gitc/gatk4/gatk; \
    \
    rm -r /usr/gitc/temp \
    ; \
# Symlink python to python3
    ln -s python3 /usr/bin/python;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
