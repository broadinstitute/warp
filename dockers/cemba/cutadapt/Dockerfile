FROM python:2

LABEL maintainer="DSP Data Centers <monster@broadinstitute.org>"

ARG CUTADAPT_VERSION

# get cutadapt
RUN pip install --upgrade pip && \
      pip install cutadapt==${CUTADAPT_VERSION}
