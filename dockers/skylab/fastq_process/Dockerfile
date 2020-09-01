FROM ubuntu:16.04
LABEL maintainer=" Kishori Jishu Xu <jishuxu@broadinstitute.org> "\
      software="fastqproces" \
      version="1.0.0" \
      description="Processes FASTQ files for 10x genomic inputs" \
      website="TBA"


# install additional python packages
#Install wget, unzip
RUN apt update && apt install -y \
  wget && \
  apt-get install -y make build-essential  zlib1g-dev 

#  liblzma-dev \
#  libbz2-dev \
#  build-essential \
#  zlib1g-dev 

RUN wget https://github.com/kishorikonwar/fastqpreprocessing/archive/1.0.0.tar.gz

RUN tar -zxvf 1.0.0.tar.gz && \
    cd fastqpreprocessing-1.0.0 &&\
    ls &&\
    make -C libStatGen &&\
    mkdir src/obj &&\
    make -C src/ 

RUN cp fastqpreprocessing-1.0.0/src/fastqprocess /usr/local/bin/


WORKDIR /
