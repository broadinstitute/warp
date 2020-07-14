FROM python:3.6.2

LABEL maintainer="Nick Barkas <nbarkas@broadinstitute.org>" \
    software="umi_tools" \
    version="0.5.5" \
    description="tools for extraction correction, deduplication and counting of UMIs" \
    website="https://github.com/CGATOxford/UMI-tools"

RUN git clone https://github.com/CGATOxford/UMI-tools.git
WORKDIR UMI-tools
RUN git checkout tags/0.5.5
RUN pip install .

RUN mkdir /root/tools
COPY getUntaggedReads /root/tools

ENV PATH="/root/tools/:$PATH"

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