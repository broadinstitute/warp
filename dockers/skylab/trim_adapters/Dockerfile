FROM ubuntu:16.04
LABEL software="ea-utils package" \
      version="1.04.807" \
      description="Command-line tools for processing biological sequencing data. Barcode demultiplexing, adapter trimming, etc. Primarily written to support an Illumina based pipeline - but should work with any FASTQs." \
      website="http://expressionanalysis.github.io/ea-utils/"

# Install compiler 
RUN apt update && apt install -y \
   zlib1g-dev \
   cmake \
   automake \
   wget \
   build-essential \	
   libgsl0-dev 

# Install ea-utils 
WORKDIR /usr/local/ 
RUN wget "https://github.com/ExpressionAnalysis/ea-utils/archive/refs/tags/1.04.807.tar.gz"
RUN tar -xzvf 1.04.807.tar.gz
WORKDIR /usr/local/ea-utils-1.04.807/clipper 
RUN make
RUN make install
ENV PATH /usr/local/ea-utils-1.04.807/clipper/:$PATH
# Cleanup
RUN apt-get clean
