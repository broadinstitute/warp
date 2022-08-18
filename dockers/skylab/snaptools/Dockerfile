FROM --platform=linux/amd64 python:3.7.2

ENV TERM=xterm-256color \
        TINI_VERSION=v0.19.0 \
        PATH=$PATH:/usr/gitc

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

RUN set -eux; \
    pip install --upgrade pip; \
    apt-get update; \
    apt-get install -y bash wget \
    ; \
    # Install snaptools
    pip install snaptools; \
    mkdir -p /home/r3fang/anaconda2/bin/; \
    ln -s `which python` /home/r3fang/anaconda2/bin/python \
    ; \
    # Install TINI
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]