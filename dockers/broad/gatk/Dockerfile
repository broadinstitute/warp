FROM adoptopenjdk/openjdk8:alpine-slim

ARG GATK4_VERSION=4.1.8.0 \
        GATK3_VERSION=3.5

ENV TERM=xterm-256color \
    GATK4_URL=https://github.com/broadinstitute/gatk/archive/refs/tags/${GATK4_VERSION}.tar.gz \
    GATK3_URL=https://broadinstitute.jfrog.io/artifactory/libs-release-local/org/broadinstitute/gatk/gatk-utils/${GATK3_VERSION}/gatk-utils-${GATK3_VERSION}.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        GATK4_VERSION=${GATK4_VERSION} \
        GATK3_VERSION=${GATK3_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            bash \
            tini \ 
            wget \
            python3 \
            py3-pip \
        ; \
        echo 'alias python="/usr/bin/python3"' >> ~/.bashrc; \
            source ~/.bashrc; \
        \
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
        tar xf ${GATK4_VERSION}.tar.gz; \
        mv gatk-${GATK4_VERSION}/gatk /usr/gitc/gatk4/; \
        \
        rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
