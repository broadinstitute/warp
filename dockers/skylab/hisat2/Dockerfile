FROM ubuntu:16.04
LABEL MAINTAINER="Jishu Xu <jishuxu@broadinstitute.org>"
LABEL software="HISAT2"
LABEL version="2-2.1.0"
LABEL description="RNA-seq aligner"
LABEL website="https://ccb.jhu.edu/software/hisat2/index.shtml"

RUN mkdir -p /opt/tools/
WORKDIR /opt/tools

RUN \
 apt update && \ 
 apt install -y \
  liblzma-dev \
  libbz2-dev \
  cmake \
  automake \
  curl \
  libboost-all-dev \
  libcurl4-openssl-dev \
  wget \
  build-essential \	
  gcc-multilib \
  zlib1g-dev \
  libxml2-dev \
  libncurses5-dev \
  zip unzip \
  git \
  r-base \
  r-base-core \
  r-base-dev
  
RUN \
  wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip && \
  unzip hisat2-2.1.0-source.zip && \
  cd hisat2-2.1.0 && \
  make && \
  cp hisat2* /usr/local/bin

COPY  ./*.sh hisat2-2.1.0/
# set ENV PATH, run hisat2 by type 'HISAT2' from cmd
ENV PATH="/opt/tools/hisat2-2.1.0:${PATH}"

# Install samtools
WORKDIR /usr/local/samtools
ADD https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 .

RUN tar -xvf samtools-1.6.tar.bz2 && \
    rm samtools-1.6.tar.bz2 && \
    cd samtools-1.6 && \
    ./configure --prefix=/usr && \
    make && \
    make install

# gffread is gtf/gtt tools 
WORKDIR /opt/tools/gffread
RUN \
  git clone https://github.com/gpertea/gclib && \
  git clone https://github.com/gpertea/gffread && \
  cd gffread && \
  make && \ 
  cp gffread /usr/local/bin/
WORKDIR /opt/tools
