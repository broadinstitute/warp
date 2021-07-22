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
  && rm -rf /var/lib/apt/lists/*

#Should be same eagle version as Michigan imputation server
ARG EAGLE_VERSION=2.4
# RUN wget https://github.com/genepi/imputationserver/raw/master/files/bin/eagle
RUN wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/old/Eagle_v${EAGLE_VERSION}.tar.gz
RUN tar xf Eagle_v2.4.tar.gz
RUN mv Eagle_v2.4/eagle .
RUN rm -r Eagle_v2.4.tar.gz
RUN ["chmod", "+x", "eagle"]

RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz
