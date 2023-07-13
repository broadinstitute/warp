FROM python:3.7.2

RUN set -eux; \
    apt-get update; \
    apt-get install -y git-all; \
    apt-get install wget; \
    apt-get install make; \
    git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n; \
    cd hisat-3n; \
    git checkout -b hisat-3n origin/hisat-3n; \
    make

