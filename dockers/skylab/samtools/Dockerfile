FROM ubuntu:16.04

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="samtools" \
  version="1.6" \
  description="processing sequence alignments in SAM and BAM formats" \
  website="https://samtools.github.io"

RUN apt update && \
    apt install -y \
      wget \
      bzip2 \
      g++ \
      cmake \
      curl \
      libncurses5-dev \
      zlib1g-dev \
      libbz2-dev \
      zip \
      unzip \
      liblzma-dev \
      openssl \
      libcurl4-openssl-dev \
      libssl-dev

WORKDIR /usr/local/samtools
ADD https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 .

RUN tar -xvf samtools-1.6.tar.bz2 && \
    rm samtools-1.6.tar.bz2 && \
    cd samtools-1.6 && \
    ./configure --prefix=/usr && \
    make && \
    make install
