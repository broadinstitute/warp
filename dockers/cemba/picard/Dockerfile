FROM alpine:3.8

LABEL maintainer="DSP DATA CENTERS <monster@broadinstitute.org>"

ARG PICARD_VERSION

# ca-certificates (for http requests)
# libc6 (for Intel GKL)
# R for picard
RUN apk --update add \
      bash \
      openjdk8 \
      ca-certificates \
      libc6-compat \
      unzip \
      curl \
      R \
      R-dev \ 
      wget && \
      update-ca-certificates

# get picard 
# older versions/releases have a zip not a jar (2.5.0) 
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar && \
      mkdir picard-tools && \
      mv picard.jar picard-tools/

# TODO: Testing with the cloud jar. If that works we can remove this from here and the build script