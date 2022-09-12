FROM --platform=linux/amd64 python:3.7.2

ENV TERM=xterm-256color \
        TINI_VERSION=v0.19.0\
        PATH=$PATH:/usr/gitc

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

# Copy scripts and dependencies into image
COPY tools .
COPY requirements.txt .

RUN set -eux; \
    pip install --upgrade pip; \
    apt-get update; \
    apt-get install -y bash wget; \
    python -m pip install git+https://github.com/HumanCellAtlas/sctools.git#egg=sctools; \
    pip3 install -r requirements.txt; \
    # Install TINI
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--"]