FROM quay.io/broadinstitute/bowtie2:2.3.4.3

LABEL maintainer="DSP Data Centers <monster@broadinstitute.org>"

ARG BISMARK_VERSION

# get bismark (wrapper for bowtie2) 
# older versions/releases have v before version (0.15.0)
RUN wget https://github.com/FelixKrueger/Bismark/archive/${BISMARK_VERSION}.tar.gz && \
      tar xf ${BISMARK_VERSION}.tar.gz && \
      rm ${BISMARK_VERSION}.tar.gz && \
      mv /Bismark-${BISMARK_VERSION} /Bismark
