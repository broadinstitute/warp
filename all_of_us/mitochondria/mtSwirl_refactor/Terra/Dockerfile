# Docker image for Step 3 (combine VCFs + homref imputation from coverage.h5)
#
# This image MUST include:
#   - Hail + Spark runtime (we deliberately inherit from the current AoU mito image)
#   - mtSwirl source code at /opt/mtSwirl so WDL can call:
#       python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/combine_vcfs_and_homref_from_covdb.py
#
# Build context should be the mtSwirl repo root:
#   docker build -f generate_mtdna_call_mt/Terra/Dockerfile -t aou-mitochondrial-combine-vcfs-covdb:dev .

FROM us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0

ENV PYTHONUNBUFFERED=1 \
    PYTHONPATH=/opt/mtSwirl

WORKDIR /opt/mtSwirl

# Add Google Cloud SDK tooling so WDL tasks can manually download gs:// inputs
# (e.g. shard MT tarballs) at runtime without relying on Cromwell localization.
# We install the minimal packages needed for `gcloud storage cp` and `gsutil cp`.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        gnupg \
        apt-transport-https \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && apt-get update \
    && apt-get install -y --no-install-recommends google-cloud-cli \
    && rm -rf /var/lib/apt/lists/*

# Step 3 needs to read the v2 coverage DB (`coverage.h5`). The base image does not
# necessarily include h5py, so install it explicitly.
#
# Notes:
# - We keep this minimal to avoid bloating the Hail image.
# - h5py provides manylinux wheels for linux/amd64.
RUN python3 -m pip install --no-cache-dir \
    h5py \
    gcsfs

# Layer in the mtSwirl code. We copy the full package to avoid missing module imports.
# NOTE: build context should be the mtSwirl_refector directory.
COPY ./ /opt/mtSwirl/generate_mtdna_call_mt/

# Keep entrypoint neutral.
ENTRYPOINT ["python3"]
