FROM ubuntu:16.04
LABEL maintainer="Kylee Degatano<kdegatano@broadinstitute.org> "\
      software="STAR" \
      version="2.7.9a" \
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
WORKDIR /usr/local/
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz && \
  tar -xf 2.7.9a.tar.gz
RUN chmod +x /usr/local/STAR-2.7.9a/bin/Linux_x86_64_static/STAR
ENV PATH /usr/local/STAR-2.7.9a/bin/Linux_x86_64_static/:$PATH
RUN cp /usr/local/STAR-2.7.9a/bin/Linux_x86_64_static/STAR /usr/local/bin

WORKDIR /
RUN rm -rf /build
