FROM marketplace.gcr.io/google/debian9
MAINTAINER DSDE <lantern@broadinstitute.org>

ARG PICARD_PUBLIC_VERSION=2.23.8
ARG GATK35_VERSION=3.5-0-g36282e4
ARG GATK4_VERSION=4.1.8.0
ARG SAMTOOLS_VER=1.11
ARG BWA_VER=0.7.15.r1140
ARG TABIX_VER=0.2.5_r1005
ARG BGZIP_VER=1.3
ARG SVTOOLKIT_VER=2.00-1650

ENV TERM=xterm-256color
ENV DOCKER_FIX='                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        '

LABEL GOTC_PICARD_VER=$PICARD_PUBLIC_VERSION
LABEL GOTC_GATK35_VER=$GATK35_VERSION
LABEL GOTC_GATK4_VER=$GATK4_VERSION
LABEL GOTC_SAMTOOLS_VER=$SAMTOOLS_VER
LABEL GOTC_BWA_VER=$BWA_VER
LABEL GOTC_TABIX_VER=$TABIX_VER
LABEL GOTC_BGZIP_VER=$BGZIP_VER
LABEL GOTC_SVTOOLKIT_VER=$SVTOOLKIT_VER

# Change working directory to /usr/gitc/
WORKDIR /usr/gitc

# Install java8, python2, python3, and R
RUN apt-get update && \
	apt-get upgrade -y && \
	apt-get install -y openjdk-8-jdk && \
	apt-get install -y python && \
	apt-get install -y python3-pip && \
	apt-get install -y r-base && \
	apt-get install wget

# Install ggplot2
RUN echo 'install.packages(c(\"ggplot2\"), repos=\"http://cran.us.r-project.org\", dependencies=TRUE)' > /tmp/packages.R && \
    Rscript /tmp/packages.R

# Install samtools
RUN wget 'https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2' && \
  tar xf samtools-${SAMTOOLS_VER}.tar.bz2 && \
  rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
  cd samtools-${SAMTOOLS_VER} && \
  ./configure && \
  make && \
  make install && \
  cd ../ && \
  rm -r samtools-${SAMTOOLS_VER} # Don't need the source directory now that we've installed the binary

RUN export SHORT_BWA_VER="${BWA_VER:0:6}" && \
  wget 'https://github.com/lh3/bwa/releases/download/v${SHORT_BWA_VER}/bwakit-${SHORT_BWA_VER}_x64-linux.tar.bz2' -O bwakit-${SHORT_BWA_VER}.tar.bz2 && \
  tar xf bwakit-${SHORT_BWA_VER}.tar.bz2 && \
  rm bwakit-${SHORT_BWA_VER}.tar.bz2 && \
  mv bwa.kit/bwa bwa && \
  rm -rf bwa.kit

# Copy everything in working dir outside container to /usr/gitc
COPY . .
