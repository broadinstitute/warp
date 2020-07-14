# Start with an ubuntu system and update it
FROM ubuntu:18.04

# Avoid interaction during package install
ENV DEBIAN_FRONTEND=noninteractive

# Image label
LABEL maintainer="Nick Barkas <nbarkas@broadinstititute.org>" \
  software="DropletUtils" \
  version="1.2.1-0.1.0" \
  description="Bioconductor DropletUtils Package with the suitable version of R" \
  website="https://bioconductor.org/packages/release/bioc/html/DropletUtils.html"

# Enable source repositories to install deps for R and update the apt-get list
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && apt-get update

# Install R dependencies, wget, Python and numpy
RUN apt-get install -y \
    build-essential \
    wget \
    python \
    python-pip \
    && apt-get build-dep -y r-base

RUN pip install numpy

# Download and unzip R Build R and add to PATH
WORKDIR /root/
RUN wget https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz \
    && tar xzf R-3.5.1.tar.gz \
    && cd "R-3.5.1" \
    && ./configure --prefix=/opt/R/3.5.1/ --enable-R-shlib --with-blas --with-lapack \
    && make \
    && make install \
    && cd /root/ \
    && rm -r R-3.5.1 \
    && rm R-3.5.1.tar.gz

ENV PATH="/opt/R/3.5.1/bin/:$PATH"

# Install R dependencies and clean up
COPY installRDeps.R /root/
RUN Rscript /root/installRDeps.R \
    && rm installRDeps.R

# Add EmptyDrop scripts
RUN mkdir /root/tools
COPY emptyDropsWrapper ./tools/emptyDropsWrapper
COPY npz2rds ./tools/npz2rds

ENV PATH="/root/tools/emptyDropsWrapper/:/root/tools/npz2rds/:$PATH"
