FROM ubuntu:16.04
LABEL maintainer=" Kishori Konwar <kkonwar@broadinstitute.org> "\
      software="fastqproces" \
      version="1.0.0" \
      description="Processes FASTQ files for 10x genomic inputs" \
      website="TBA"

# install additional python packages
#Install wget, unzip
RUN apt update && apt install -y \
  wget zip && \
  apt-get install -y make build-essential  zlib1g-dev 

#  liblzma-dev \
#  libbz2-dev \
#  build-essential \
#  zlib1g-dev 

RUN wget https://github.com/HumanCellAtlas/sctools/archive/kmk-fastqprocessing.zip 

RUN unzip kmk-fastqprocessing.zip && \
    cd sctools-kmk-fastqprocessing/fastqpreprocessing &&\
    make -C libStatGen &&\
    mkdir src/obj &&\
    make -C src/ 

RUN cp sctools-kmk-fastqprocessing/fastqpreprocessing/src/fastqprocess /usr/local/bin/


WORKDIR /
