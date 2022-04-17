# Minimac image uses bcftools
FROM us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.5-1.10.2-0.1.16-1649948623

ARG MINIMAC4_VERSION=1.0.2

ENV TERM=xterm-256color

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        MINIMAC4_VERSION=${MINIMAC4_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            bash \
    ; \
# Install Minimac4
    wget https://github.com/statgen/Minimac4/releases/download/v${MINIMAC4_VERSION}/minimac4-${MINIMAC4_VERSION}-Linux.sh; \
    yes | bash minimac4-${MINIMAC4_VERSION}-Linux.sh --prefix=/usr/temp; \
    \
    mv /usr/temp/minimac4-${MINIMAC4_VERSION}-Linux/bin/minimac4 /usr/gitc; \
    rm minimac4-${MINIMAC4_VERSION}-Linux.sh

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
