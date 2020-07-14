FROM openjdk:8-jre

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="dropseqtools" \
  version="1.12" \
  description="tools for manipulation of drop-seq data and BAM files" \
  website="http://mccarrolllab.com/dropseq/"

RUN apt update && apt install -y \
  curl \
  unzip

RUN  apt install -y python

RUN curl http://mccarrolllab.com/download/922/ >> Drop-seq_tools-1.12.zip && \
  unzip Drop-seq_tools-1.12.zip && \
  cp -r Drop-seq_tools-1.12/* /usr/local/bin/
