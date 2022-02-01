FROM adoptopenjdk/openjdk8:alpine-slim

ARG FGBIO_VERSION=1.4.0

ENV TERM=xterm-256color \
        FGBIO_URL=https://github.com/fulcrumgenomics/fgbio/releases/download/${FGBIO_VERSION}/fgbio-${FGBIO_VERSION}.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        FGBIO_VERSION=${FGBIO_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk --no-cache upgrade; \
        apk add --no-cache \ 
            bash \
            wget \
            tini \
    ; \
# Install fgbio jar
    wget ${FGBIO_URL} -O fgbio.jar 

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]