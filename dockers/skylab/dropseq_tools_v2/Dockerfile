FROM openjdk:8-jre

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>" \
  software="dropseqtools" \
  version="2.3.0" \
  description="tools for manipulation of drop-seq data and BAM files" \
  website="http://mccarrolllab.com/dropseq/"

RUN apt update && apt install -y \
  curl \
  unzip

RUN  apt install -y python

RUN curl -L https://github.com/broadinstitute/Drop-seq/releases/download/v2.3.0/Drop-seq_tools-2.3.0.zip > Drop-seq_tools-2.3.0.zip && \
  unzip Drop-seq_tools-2.3.0.zip && \
  chmod +x Drop-seq_tools-2.3.0/* && \
  cp -r Drop-seq_tools-2.3.0/* /usr/local/bin/ && \
  rm Drop-seq_tools-2.3.0.zip && \
  rm -rf Drop-seq_tools-2.3.0/
