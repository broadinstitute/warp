FROM --platform=linux/amd64 ubuntu:16.04

ARG BWA_VERSION=0.7.17

ENV TERM=xterm-256color \
        TINI_VERSION=v0.19.0 \
        BWA_URL=https://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2/download

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        BWA_VERSION=${BWA_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
        apt update; \
        apt install -y \
            bzip2 \
            zlib1g \
            zlib1g-dev \
            wget \
    ; \
# Install BWA \
    mkdir temp; \
    cd temp; \
    wget -O "bwa.tar.bz2" ${BWA_URL}; #\
    tar -xvjf bwa.tar.bz2; \
    cd bwa-${BWA_VERSION}; \
    make; \
    chmod +x ./bwa; \
    cp ./bwa /usr/local/bin; \
    rm -r /usr/gitc/temp \
#    ; \
# Install TINI
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini \
    ; \

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--"]

# Install system requirements
#RUN apt update  && \
#    apt install -y git python python-pip zlib1g zlib1g-dev wget

# Install SnapTools
#RUN mkdir /install && \
#    cd /install && \
#    git clone https://github.com/r3fang/SnapTools.git && \
#    cd SnapTools && \
#    pip install .

#RUN cd /install && \
#    wget -O "bwa-0.7.17.tar.bz2" "https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download" && \
#    tar xvjf bwa-0.7.17.tar.bz2 && \
#    cd bwa-0.7.17 && \
#    make && \
#    mkdir /tools/ && \
#    cp bwa /tools/

# This fixes an error in the code of snaptools that
# references an absolute path for the python binary
#RUN mkdir -p /home/r3fang/anaconda2/bin/ && \
#    ln -s `which python` /home/r3fang/anaconda2/bin/python

