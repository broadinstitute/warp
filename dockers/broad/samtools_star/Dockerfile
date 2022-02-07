FROM us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616

ARG STAR_VERSION=2.7.10a

ENV TERM=xterm-256color \
        STAR_URL=https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        STAR_VERSION=${STAR_VERSION}

WORKDIR /usr/gitc

RUN set -eux \
    ; \
# Install STAR
    mkdir temp; \
    cd temp; \
    \
    wget ${STAR_URL}; \
    tar -xzf ${STAR_VERSION}.tar.gz; \
    cd STAR-${STAR_VERSION}/bin/Linux_x86_64_static; \
    \
    chmod +x ./STAR; \
    cp ./STAR /usr/local/bin; \
    \
    rm -r /usr/gitc/temp

# Set tini as default entrypoint 
ENTRYPOINT [ "/sbin/tini", "--" ]
