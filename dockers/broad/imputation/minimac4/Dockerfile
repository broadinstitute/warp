# Minimac image uses bcftools
FROM us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.1-1.10.2-0.1.16-1637264085

ARG MINIMAC4_VERSION=1.0.2

ENV TERM=xterm-256color

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        MINIMAC4_VERSION=${MINIMAC4_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            bash \
            cmake \
    ; \
# Install python dependencies
    pip3 install cget \
    ; \
# Install Minimac4
    wget https://github.com/statgen/Minimac4/archive/v${MINIMAC4_VERSION}.tar.gz; \
    tar xf v${MINIMAC4_VERSION}.tar.gz; \
    \
    cd Minimac4-${MINIMAC4_VERSION}; \
    bash install.sh; \
    cp release-build/minimac4 ..; \
    cd ..;\
    \
    rm -r Minimac4-${MINIMAC4_VERSION}; \
    rm v${MINIMAC4_VERSION}.tar.gz

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
