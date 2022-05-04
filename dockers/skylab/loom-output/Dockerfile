FROM python:3.7.2

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>"

RUN pip install --upgrade pip

RUN apt-get update && apt-get install wget

RUN python -m pip install git+https://github.com/HumanCellAtlas/sctools.git#egg=sctools

COPY requirements.txt .
RUN pip3 install -r requirements.txt

RUN mkdir /tools
WORKDIR /tools

COPY create_loom_optimus.py .
COPY create_loom_ss2.py .
COPY loomCompare.py .
COPY ss2_loom_merge.py .
COPY create_snss2_counts_csv.py .
COPY create_loom_snss2.py .
COPY create_snrna_optimus.py .
COPY create_snrna_optimus_counts.py .
