# Adding a platform tag to ensure that images built on ARM-based machines (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
FROM --platform="linux/amd64" ubuntu:18.04

# Avoid interaction during package install
ENV DEBIAN_FRONTEND=noninteractive

# Image label
LABEL maintainer="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
  software="DropletUtils" \
  version="1.2.1-0.1.0" \
  description="Bioconductor DropletUtils Package with the suitable version of R" \
  website="https://bioconductor.org/packages/release/bioc/html/DropletUtils.html"

# Download and unzip R Build R and add to PATH
WORKDIR /usr/local/
COPY installRDeps.R /usr/local/
# Enable source repositories to install deps for R and update the apt-get list
# Install R dependencies, wget, Python and numpy. Cleanup afterwards.
# NOTE: && is used to chain commands so it fails fast.
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && apt-get update \
    && apt-get install -y \
    build-essential \
    wget \
    python \
    python-pip \
    && apt-get build-dep -y r-base \
    && pip install numpy \
    && wget https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz \
    && tar xzf R-3.5.1.tar.gz \
    && cd "R-3.5.1" \
    && ./configure --prefix=/opt/R/3.5.1/ --enable-R-shlib --with-blas --with-lapack \
    && make \
    && make install \
    && cd /usr/local/ \
    && rm -r R-3.5.1 \
    && rm R-3.5.1.tar.gz \
    && /opt/R/3.5.1/bin/Rscript /usr/local/installRDeps.R \
    && rm installRDeps.R \
    && mkdir /usr/local/tools

# Add R installation directory to path
ENV PATH="/opt/R/3.5.1/bin/:$PATH"

# Add EmptyDrop scripts
COPY emptyDropsWrapper ./tools/emptyDropsWrapper
COPY npz2rds ./tools/npz2rds

ENV PATH="/usr/local/tools/emptyDropsWrapper/://usr/local/tools/npz2rds/:$PATH"
