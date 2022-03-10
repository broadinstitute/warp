FROM python:3.7.2

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>"

RUN pip install --upgrade pip

RUN mkdir /tools
WORKDIR /tools

COPY convert-vcf-to-allc.py .
