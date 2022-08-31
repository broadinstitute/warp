FROM --platform=linux/amd64 python:3.6.2-alpine

ARG SUBREAD_VERSION="2.0.1"

ENV TERM=xterm-256color \
        SUBREAD_URL="https://downloads.sourceforge.net/project/subread/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-source.tar.gz" \
        PATH=$PATH:/usr/gitc/subread-${SUBREAD_VERSION}-source/bin

LABEL MAINTANER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
      SUBREAD_VERSION=${SUBREAD_VERSION}

WORKDIR /usr/gitc

COPY requirements.txt .
COPY remove-reads-on-junctions.py .

RUN set -eux; \
        apk add --no-cache \
            bash \
            build-base \
            bzip2-dev \
            tini \
            wget \
            xz-dev \
            zlib-dev; \
    python3 -m pip install --upgrade pip; \
    pip3 install -r requirements.txt \
    ; \
# Install subread
    wget --no-check-certificate ${SUBREAD_URL}; \
    tar -xzvf subread-${SUBREAD_VERSION}-source.tar.gz; \
    cd /usr/gitc/subread-${SUBREAD_VERSION}-source/src; \
    make -f Makefile.Linux;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]