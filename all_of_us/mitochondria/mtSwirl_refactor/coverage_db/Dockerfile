# Minimal image to run the v2 coverage DB builder (HDF5 + exact summary metrics)
#
# Build context should be the repo root so we can COPY generate_mtdna_call_mt/.
# Example:
#   docker build -f generate_mtdna_call_mt/coverage_db/Dockerfile -t mt-coverage-db:dev .

FROM python:3.12-slim

# Keep Python output unbuffered for better Cromwell logs
ENV PYTHONUNBUFFERED=1 \
  PIP_NO_CACHE_DIR=1 \
  PYTHONPATH=/opt/mtSwirl

WORKDIR /opt/mtSwirl

# System deps:
# - libhdf5 via manylinux wheels works without apt-get, but having a compiler toolchain
#   can help if pip ever has to build from source.
# - gzip/tar/coreutils are already available in slim, but we install a minimal set.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
       ca-certificates \
  && rm -rf /var/lib/apt/lists/*

# Python deps. Pin loosely to avoid unnecessary churn but keep reproducibility.
# If you need stricter pinning for production, we can add exact versions.
RUN python -m pip install --upgrade pip \
  && python -m pip install \
       numpy \
       pandas \
  h5py \
  gcsfs

# Copy only the Python package code needed.
COPY ./ /opt/mtSwirl/generate_mtdna_call_mt/

# Default entrypoint: keep it neutral so WDL can invoke any module.
ENTRYPOINT ["python"]
