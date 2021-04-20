FROM ubuntu:16.04
LABEL maintainer=" Jishu Xu <jishuxu@broadinstitute.org> "\
      software="STAR" \
      version="2.5.3a" \
      description="RNA-seq aligner" \
      website="https://github.com/alexdobin/STAR"
RUN mkdir build 
WORKDIR build 
# install additional python packages
#Install wget, unzip
RUN apt update && apt install -y \
  liblzma-dev \
  libbz2-dev \
  cmake automake \
  curl \
  libboost-all-dev \
  wget \
  build-essential \
  gcc-multilib \
  zlib1g-dev \
  libxml2-dev \
  libncurses5-dev \
  r-base \
  r-base-core \
  r-base-dev 

RUN wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
  tar -xf 2.5.3a.tar.gz
WORKDIR STAR-2.5.3a/bin/Linux_x86_64_static
RUN cp STAR /usr/local/bin

WORKDIR /
RUN rm -rf /build
