#!/bin/bash

TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

set -e

TMPDIR=$(mktemp -t gotc-dockerXXXXXX -d)

DOCKER_VERSION="2.4.7"
DOCKER_IMAGE_TAG="$DOCKER_VERSION-$TIMESTAMP"

PICARD_PRIVATE_VERSION="1.1448"
PICARD_PUBLIC_VERSION="2.23.8"
GATK35_VERSION="3.5-0-g36282e4"
GATK4_VERSION="4.1.8.0"
SAMTOOLS_VER="1.11"
BWA_VER="0.7.15.r1140"
SHORT_BWA_VER="${BWA_VER:0:6}"
K8_VER="null"
BWA_POSTALT_SCRIPT_VER="null"
TABIX_VER="0.2.5_r1005"
BGZIP_VER="1.3"
SVTOOLKIT_VER="2.00-1650"

PICARD="/seq/software/picard-public/${PICARD_PUBLIC_VERSION}/picard.jar"
GATK35="/seq/software/gotc/gatk/GenomeAnalysisTK-${GATK35_VERSION}/GenomeAnalysisTK-${GATK35_VERSION}.jar"
GATK4="/seq/software/gotc/gatk/gatk4/gatk-${GATK4_VERSION}/"
TABIX="/seq/software/picard/${PICARD_PRIVATE_VERSION}/3rd_party/tabix/tabix"
BGZIP="/seq/software/gotc/3rd_party/bgzip/bgzip"
SVTOOLKIT="/seq/software/gotc/svtoolkit/svtoolkit2.00/"

scp -T vpicard05:"$PICARD $GATK35 $TABIX $BGZIP" "${TMPDIR}"/
scp -r vpicard05:"$SVTOOLKIT" ${TMPDIR}/
scp -r vpicard05:"$GATK4" ${TMPDIR}/gatk4/

cd ${TMPDIR}

mv GenomeAnalysisTK-${GATK35_VERSION}.jar GATK35.jar

echo "FROM marketplace.gcr.io/google/debian9
MAINTAINER DSDE <lantern@broadinstitute.org>

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

RUN wget 'https://github.com/lh3/bwa/releases/download/v${SHORT_BWA_VER}/bwakit-${SHORT_BWA_VER}_x64-linux.tar.bz2' -O bwakit-${SHORT_BWA_VER}.tar.bz2 && \
  tar xf bwakit-${SHORT_BWA_VER}.tar.bz2 && \
  rm bwakit-${SHORT_BWA_VER}.tar.bz2 && \
  mv bwa.kit/bwa bwa && \
  rm -rf bwa.kit

# Copy everything in working dir outside container to /usr/gitc
COPY . .
" > Dockerfile

echo -e "$DOCKER_IMAGE_TAG\t$PICARD_PUBLIC_VERSION\tn/a\t$GATK35_VERSION\tn/a\t$GATK4_VERSION\t$SAMTOOLS_VER\t$BWA_VER\t$K8_VER\t$BWA_POSTALT_SCRIPT_VER\t$TABIX_VER\t$BGZIP_VER\t$SVTOOLKIT_VER" >> ${DIR}/build_docker_version.tsv

# Tagged with the picard release version the jars and binaries were taken from.
docker build -t us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:${DOCKER_IMAGE_TAG} .
docker push us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:${DOCKER_IMAGE_TAG}

rm -rf ${TMPDIR}
