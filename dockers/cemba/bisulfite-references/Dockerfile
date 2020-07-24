FROM continuumio/miniconda3:4.5.11

LABEL maintainer="DSP Data Centers <monster@broadinstitute.org>"

ARG SCRIPT_DIR

# use the python script to convert c's to t's and g's to a's
COPY ${SCRIPT_DIR}/build_bisulfite_references.py .
