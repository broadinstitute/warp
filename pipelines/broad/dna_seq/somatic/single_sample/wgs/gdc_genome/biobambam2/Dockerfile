FROM ubuntu:bionic-20180426

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

RUN set -ex \
    && apt-get update && apt-get install -y --no-install-recommends \
        make \
        curl \
        g++ \
        zlib1g-dev \
        ca-certificates \
        pkg-config \
        libgomp1 \
        
    && curl https://gitlab.com/german.tischler/libmaus2/repository/archive.tar.gz\?ref\=2.0.483-release-20180507173657 -o libmaus2.tar.gz \
    && tar -xzvf libmaus2.tar.gz \
    && cd libmaus2-2.0.483-release-20180507173657* \
    && ./configure \
    && make -j4 \
    && make install \
    && cd .. \
    && rm -rf libmaus2* \
    && curl https://gitlab.com/german.tischler/biobambam2/repository/archive.tar.gz\?ref\=2.0.87-release-20180301132713 -o biobambam2.tar.gz \
    && tar -xzvf biobambam2.tar.gz \
    && cd biobambam2-2.0.87-release-20180301132713* \
    && export LIBMAUSPREFIX=/usr/local \
    && ./configure --with-libmaus2=${LIBMAUSPREFIX} \
    && make -j4 \
    && make install \
    && cd .. \
    && rm -rf biobambam2* \
    && apt-get remove --purge -y \
        make \
        curl \
        ca-certificates \
        pkg-config \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*