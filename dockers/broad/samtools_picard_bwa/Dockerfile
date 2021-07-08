FROM us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616

ARG BWA_VERSION=0.7.15 \
        PICARD_PUBLIC_VERSION=2.23.8

ENV TERM=xterm-256color \
    BWA_URL=https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwakit-${BWA_VERSION}_x64-linux.tar.bz2 \
    PICARD_URL=https://github.com/broadinstitute/picard/releases/download/${PICARD_PUBLIC_VERSION}/picard.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        BWA_VERSION=${BWA_VERSION} \
        PICARD_PUBLIC_VERSION=${PICARD_PUBLIC_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            libc6-compat \
            gcompat \
        ; \
# Install BWA 
        wget ${BWA_URL} -O bwakit-${BWA_VERSION}.tar.bz2; \
        tar xf bwakit-${BWA_VERSION}.tar.bz2;  \
        rm bwakit-${BWA_VERSION}.tar.bz2; \
        mv bwa.kit/bwa bwa; \
        rm -rf bwa.kit \
        ; \
# Install Picard
        wget ${PICARD_URL}

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--"