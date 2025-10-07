FROM python:3.7.7

LABEL maintainer="kylee@broadinstitute.org" \
  software="python 3.7.7" \
  description="python 3.7.7 scripts to manage terra workflow with fiss"


RUN apt-get update && apt-get install -y tar && apt-get install -y vim

RUN mkdir /terra-scripts/

COPY . /terra-scripts 

RUN pip3 install pandas

RUN cd /terra-scripts &&\
    wget https://github.com/broadinstitute/fiss/archive/v0.16.30.tar.gz &&\
    tar -zxvf v0.16.30.tar.gz &&\
    cd fiss-0.16.30 &&\
    pip3 install . 

WORKDIR /terra-scripts


