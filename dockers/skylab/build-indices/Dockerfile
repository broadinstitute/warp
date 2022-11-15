# Adding a platform tag to ensure that images built on ARM-based machines (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
# Also upgraded to a more recent Debian version (Bullseye) so tini is part of apt.
FROM --platform="linux/amd64" python:3.6-bullseye

# Add metadata labels
LABEL maintainer="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
      software="subread package" \
      version="2.0.2" \
      description="RNA-seq high-performance read alignment, quantification and mutation discovery" \
      website="http://subread.sourceforge.net/"

# Set parameters
ARG STAR_VERSION=2.7.10a

# Go into the directory where STAR will be installed
WORKDIR /usr/local/
# Copy the requirements file for pip3
COPY requirements.txt .
# Install compiler, subread, and other prerequisite libraries
# NOTE: && is used so it fails fast.
RUN apt-get update -y --fix-missing && apt-get install -y \
  wget \
  tini && \
  pip3 install --upgrade pip && \
  pip3 install -r requirements.txt && \
  wget https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz && \ 
  tar -xf ${STAR_VERSION}.tar.gz && \
  chmod +x /usr/local/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR && \
  mkdir /script

# Add STAR installation directory to PATH
ENV PATH /usr/local/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/:$PATH

# Copy the script that removes alignments spanning intron-exon junctions
WORKDIR /script
COPY add-introns-to-gtf.py .
COPY modify_gtf.py .
# Add all scripts to PATH
ENV PATH /script/:$PATH

# Set tini as default entrypoint
ENTRYPOINT ["/usr/bin/tini", "--"]