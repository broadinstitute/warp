FROM us.gcr.io/broad-dsp-gcr-public/base/jre:8-debian

ENV TERM=xterm-256color

LABEL MAINTAINER="DSDE <lantern@broadinstitute.org>" \
  GOTC_PICARD_VER=$PICARD_PUBLIC_VERSION \
  GOTC_GATK35_VER=$GATK35_VERSION \
  GOTC_GATK4_VER=$GATK4_VERSION \
  GOTC_SAMTOOLS_VER=$SAMTOOLS_VER \
  GOTC_BWA_VER=$BWA_VER \
  GOTC_TABIX_VER=$TABIX_VER \
  GOTC_BGZIP_VER=$BGZIP_VER \
  GOTC_SVTOOLKIT_VER=$SVTOOLKIT_VER

ARG PICARD_PUBLIC_VERSION=2.23.8 \
  GATK35_VERSION=3.5-0-g36282e4 \
  GATK4_VERSION=4.1.8.0 \
  SAMTOOLS_VER=1.11 \
  BWA_VER=0.7.15 \
  TABIX_VER=0.2.5_r1005 \
  BGZIP_VER=1.3 \
  SVTOOLKIT_VER=2.00-1650

# Change working directory to /usr/gitc/
WORKDIR /usr/gitc

# Install java8, python2, python3, and R
RUN apt-get update && apt-get install -y python python3-pip r-base wget && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install ggplot2
RUN echo 'install.packages("ggplot2", repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R && \
    Rscript /tmp/packages.R

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
  tar xf samtools-${SAMTOOLS_VER}.tar.bz2 && \
  rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
  cd samtools-${SAMTOOLS_VER} && \
  ./configure && \
  make && \
  make install && \
  cd ../ && \
  rm -r samtools-${SAMTOOLS_VER}

RUN wget https://github.com/lh3/bwa/releases/download/v${BWA_VER}/bwakit-${BWA_VER}_x64-linux.tar.bz2 -O bwakit-${BWA_VER}.tar.bz2 && \
  tar xf bwakit-${BWA_VER}.tar.bz2 && \
  rm bwakit-${BWA_VER}.tar.bz2 && \
  mv bwa.kit/bwa bwa && \
  rm -rf bwa.kit

# Copy everything in working dir outside container to /usr/gitc
COPY . .
