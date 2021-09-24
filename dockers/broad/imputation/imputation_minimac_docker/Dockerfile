# default ubuntu version: 20.04
ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION}

# the following argument is needed for pkg-config to
# install without hanging on user input
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y wget \
  && apt-get install -y libgomp1 \
  && apt-get install -y build-essential \
  && apt-get install -y libz-dev \
  && apt-get install -y libbz2-dev \
  && apt-get install -y liblzma-dev  \
  && apt-get install -y pkg-config \
  && apt-get install -y tabix \
  && apt-get install -y python3 \
  && apt-get install -y python3-pip \
  && apt-get install -y software-properties-common cmake \
  && rm -rf /var/lib/apt/lists/*

RUN pip3 install cget
RUN wget https://github.com/statgen/Minimac4/archive/v1.0.2.tar.gz \
  && tar xf v1.0.2.tar.gz \
  && cd Minimac4-1.0.2 \
  && bash install.sh \
  && cp release-build/minimac4 .. \
  && cd .. \
  && rm -r v1.0.2.tar.gz

RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
  && tar xf bcftools-1.10.2.tar.bz2 \
  && cd bcftools-1.10.2 \
  && ./configure \
  && make \
  && make install \
  && cd ../ \
  && rm -r bcftools-1.10.2.tar.bz2
