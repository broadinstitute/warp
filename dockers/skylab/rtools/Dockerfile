FROM --platform=linux/amd64 python:3.8-buster

ENV TERM=xterm-256color \
        PATH="$PATH:/opt/R/3.5.1/bin/:/usr/gitc/:/usr/gitc/emptyDropsWrapper/:/usr/gitc/npz2rds/" \
        R_URL="https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz"

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org"

WORKDIR /usr/gitc

COPY checkMatrix.R .
COPY installRDeps.R .
COPY emptyDropsWrapper .
COPY npz2rds .

RUN set -eux; \
        apt-get update; \
        apt-get install -y \
            build-essential \
            r-base \
            tini \
            wget; \
    chmod -R +x /usr/gitc \
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
    rm installRDeps.R;

# Set tini as default entrypoint
ENTRYPOINT ["/usr/bin/tini", "--"]