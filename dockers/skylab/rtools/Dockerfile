FROM --platform=linux/amd64 ubuntu:18.04

ENV TERM=xterm-256color \
        DEBIAN_FRONTEND=noninteractive \
        PATH="$PATH:/opt/R/3.5.1/bin/:/usr/gitc/:/usr/gitc/emptyDropsWrapper/:/usr/gitc/npz2rds/" \
        R_URL="https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz" \
        TINI_VERSION=v0.19.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

COPY checkMatrix.R .
COPY installRDeps.R .
COPY emptyDropsWrapper .
COPY npz2rds .

RUN set -eux; \
        sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list; \
        apt-get update; \
        apt-get install -y \
            build-essential \
            python \
            python-pip \
            wget; \
        apt-get build-dep -y r-base \
    ; \
    pip install numpy \
    ; \
    wget ${R_URL}; \
    tar xzf R-3.5.1.tar.gz; \
    cd "R-3.5.1"; \
    ./configure --prefix=/opt/R/3.5.1/ --enable-R-shlib --with-blas --with-lapack; \
    make; \
    make install; \
    cd /usr/gitc; \
    rm -r R-3.5.1; \
    rm R-3.5.1.tar.gz; \
    Rscript installRDeps.R; \
    rm installRDeps.R \
    ; \
# Install TINI
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]