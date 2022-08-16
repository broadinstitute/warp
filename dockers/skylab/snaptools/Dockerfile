FROM ubuntu:18.04

# Install system requirements
RUN apt update  && \
    apt install -y git python python-pip zlib1g zlib1g-dev wget

# Install SnapTools
RUN mkdir /install && \
    cd /install && \
    git clone https://github.com/r3fang/SnapTools.git && \
    cd SnapTools && \
    pip install .

RUN cd /install && \
    wget -O "bwa-0.7.17.tar.bz2" "https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download" && \
    tar xvjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    mkdir /tools/ && \
    cp bwa /tools/

# This fixes an error in the code of snaptools that
# references an absolute path for the python binary
RUN mkdir -p /home/r3fang/anaconda2/bin/ && \
    ln -s `which python` /home/r3fang/anaconda2/bin/python

ENV PATH /tools/:$PATH
