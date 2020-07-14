FROM python:3.7.2

LABEL maintainer="Mint Team <mintteam@broadinstitute.org>" \
      software="python 3.6.2" \
      description="python 3.6.2 with pysam, sctools, requests, and a basic science stack"

RUN pip3 install \
    Click==7.0 \
    numpy==1.16.2 \
    pysam==0.15.2 \
    biopython==1.73

## Install software
RUN apt-get update && \
    apt-get install -y lsb-release && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" |  tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg |  apt-key add - && \
    apt-get update && \
    apt-get install -y google-cloud-sdk

## Install latest samtools from source
RUN mkdir /tools && \
    cd /tools && \
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar xvjf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure && \
    make -j 4 && \
    cp samtools ..

## Append tools to path
ENV PATH=/tools/:${PATH}

## Copy Scripts
COPY filterFastqByReadName.py /tools/

