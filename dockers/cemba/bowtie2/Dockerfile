FROM quay.io/broadinstitute/samtools:1.9

LABEL maintainer="DSP Data Centers <monster@broadinstitute.org>"

ARG BOWTIE2_VERSION

# get bowtie2 using bioconda
RUN conda install -q -y -c bioconda \
      bowtie2==${BOWTIE2_VERSION} && \
      conda clean --all -y
