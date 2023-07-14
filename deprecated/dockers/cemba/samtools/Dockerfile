FROM continuumio/miniconda3:4.5.11

LABEL maintainer="DSP Data Centers <monster@broadinstitute.org>"

ARG SAMTOOLS_VERSION

# add channels in correct order to avoid missing libcrypto dependency
# https://github.com/bioconda/bioconda-recipes/issues/12100
RUN conda config --add channels defaults && \
      conda config --add channels bioconda && \
      conda config --add channels conda-forge && \
      conda install -q -y samtools==${SAMTOOLS_VERSION} && \
      conda clean --all -y

