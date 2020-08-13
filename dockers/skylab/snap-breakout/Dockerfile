FROM python:3.6.2

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>" \
      software="python for breakout snap step" \
      description="python for exporting snap files into csv"

RUN pip install \
    pandas==0.20.3 \
    h5py==2.9.0

RUN mkdir /tools/
COPY breakoutSnap.py /tools/
ENV PATH=/tools/:$PATH
