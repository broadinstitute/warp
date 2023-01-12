# Minimac image uses bcftools
FROM --platform=linux/amd64 us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.6-1.10.2-0.1.16-1663946207

ARG MINIMAC4_VERSION=1.0.2

ENV TERM=xterm-256color

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        MINIMAC4_VERSION=${MINIMAC4_VERSION}

WORKDIR /usr/gitc

# Install Minimac4
RUN set -eux; \
        apk upgrade; \
    wget https://github.com/statgen/Minimac4/releases/download/v${MINIMAC4_VERSION}/minimac4-${MINIMAC4_VERSION}-Linux.sh; \
    yes | bash minimac4-${MINIMAC4_VERSION}-Linux.sh --prefix=/usr/temp; \
    \
    mv /usr/temp/minimac4-${MINIMAC4_VERSION}-Linux/bin/minimac4 /usr/gitc; \
    rm minimac4-${MINIMAC4_VERSION}-Linux.sh

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--" ]
