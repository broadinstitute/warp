FROM python:3.7.2

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>"

RUN pip install --upgrade pip

RUN apt-get update

COPY requirements.txt .
RUN pip3 install -r requirements.txt

RUN mkdir /tools
WORKDIR /tools

COPY create-npz-output.py .
COPY create-merged-npz-output.py .
