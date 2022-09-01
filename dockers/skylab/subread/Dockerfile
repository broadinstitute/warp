FROM --platform=linux/amd64 python:3.6.2

ARG SUBREAD_VERSION="2.0.1"

ENV TERM=xterm-256color \
        SUBREAD_URL="https://downloads.sourceforge.net/project/subread/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-source.tar.gz" \
        TINI_VERSION=v0.19.0 \
        PATH=$PATH:/usr/gitc/subread-${SUBREAD_VERSION}-source/bin

LABEL MAINTANER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
      SUBREAD_VERSION=${SUBREAD_VERSION}

WORKDIR /usr/gitc

COPY requirements.txt .
COPY remove-reads-on-junctions.py .

RUN set -eux; \
        apt-get update --fix-missing; \
        apt-get install -y wget; \
        python3 -m pip install --upgrade pip; \
        pip3 install -r requirements.txt \
    ; \
# Install subread
    wget --no-check-certificate ${SUBREAD_URL}; \
    tar -xzvf subread-${SUBREAD_VERSION}-source.tar.gz; \
    cd /usr/gitc/subread-${SUBREAD_VERSION}-source/src; \
    make -f Makefile.Linux; \
# Install TINI \
    cd /usr/gitc; \
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]