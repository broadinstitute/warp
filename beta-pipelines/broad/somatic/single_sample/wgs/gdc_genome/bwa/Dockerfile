FROM ubuntu:bionic-20180426

RUN set -ex \
    && apt-get update && apt-get install -y --no-install-recommends \
        make \
        gcc \
        bzip2 \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        git \
        ca-certificates \
        wget \
    && mkdir bwa \
    && cd bwa \
    && git init \
    && git remote add origin https://github.com/lh3/bwa.git \
    && git fetch origin f1238714517956e32e0c9ba7eb3a882b0864bd31 \
    && git reset --hard FETCH_HEAD \
    && make \
    && mv bwa /usr/local/bin/ \
    && cd .. \
    && rm -r bwa \
    && wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 \
    && tar -xjvf samtools-1.8.tar.bz2 \
    && cd samtools-1.8 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -r samtools-1.8 \
    && apt-get remove --purge -y git ca-certificates wget $(apt-mark showauto) \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*