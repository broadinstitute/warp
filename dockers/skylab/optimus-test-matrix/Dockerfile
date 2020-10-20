# Start with an ubuntu system and update it
FROM ubuntu:18.04

# Image label
LABEL maintainer="Lantern Team <lantern@broadinstititute.org>" \
  software="Optimus Matrix Tester" 

ENV DEBIAN_FRONTEND=noninteractive

# Enable source repositories to install deps for R and update the apt-get list
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && apt-get update

# Install R dependencies, wget, Python and numpy
RUN apt-get install -y \
    build-essential \
    wget \
    python3 \
    python3-pip \
    && apt-get build-dep -y r-base

RUN pip3 install numpy

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

# Set the path
ENV PATH="/opt/R/3.5.1/bin/:$PATH"

# Install R dependencies and clean up
COPY installRDeps.R /root/
RUN Rscript /root/installRDeps.R \
    && rm installRDeps.R

# Add npz2rds tool
RUN mkdir /root/tools
COPY npz2rds ./tools/npz2rds

# Set the path
ENV PATH="/root/tools/npz2rds/:$PATH"

# Copy 
COPY checkMatrix.R /root/tools/checkMatrix.R

# Set the path
ENV PATH="/root/tools/:$PATH"