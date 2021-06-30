FROM openjdk:8-jdk-alpine

ARG GATK4_VERSION=4.1.8.0 \
        GATK35_VERSION=3.5

ENV TERM=xterm-256color \
    GATK4_URL=https://github.com/broadinstitute/gatk/archive/refs/tags/${GATK4_VERSION}.tar.gz \
    GATK35_URL=https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/gatk/gatk-utils/${GATK35_VERSION}/gatk-utils-${GATK35_VERSION}.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        GATK4_VERSION=${GATK4_VERSION} \
        GATK35_VERSION=${GATK35_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            tini \ 
            wget \
        ; \
# Install GATK
        \
        # gatk35
        wget ${GATK35_URL}; \
        mv gatk-utils-${GATK35_VERSION}.jar GATK35.jar; \
        \
        # gatk4
        wget ${GATK4_URL}; \
        tar xf ${GATK4_VERSION}.tar.gz; \
        mv gatk-${GATK4_VERSION} gatk4;


# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
