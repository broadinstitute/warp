FROM --platform=linux/amd64 python:3.7.2

ENV TERM=xterm-256color \
        PATH=$PATH:/usr/gitc

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /tools

# Copy scripts and dependencies into image
COPY tools .
COPY requirements.txt .

RUN set -eux; \
    pip install --upgrade pip; \
    apt-get update; \
    apt-get install -y bash wget; \
    python -m pip install git+https://github.com/HumanCellAtlas/sctools.git#egg=sctools; \
    pip3 install -r requirements.txt

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--"]