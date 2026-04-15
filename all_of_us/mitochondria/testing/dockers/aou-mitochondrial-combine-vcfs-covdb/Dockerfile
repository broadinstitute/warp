# Local test image for build_vcf_shard_mt and finalize_mt_with_covdb (extends production image).
#
# Changes from production:
#   1. build_vcf_shard_mt.py → wrapper that pre-localizes gs://fake-bucket/* VCFs
#      from fake-gcs-server to local disk before Hail reads them.
#   2. finalize_mt_with_covdb.py → wrapper that pre-localizes gs://fake-bucket/* BED files
#      from fake-gcs-server to local disk before Hail reads them.
#   3. /usr/local/bin/gcloud → stub that redirects `storage cp` and
#      `storage objects describe` to fake-gcs-server via curl.
#
# Build (from repo root):
#   docker build -f all_of_us/mitochondria/testing/dockers/aou-mitochondrial-combine-vcfs-covdb/Dockerfile \
#     -t aou-mitochondrial-combine-vcfs-covdb:local-test all_of_us/mitochondria/testing/dockers/aou-mitochondrial-combine-vcfs-covdb/
#
# Used by run_test.sh: the test image is shadow-tagged as the production tag
# before running miniwdl so miniwdl uses it transparently.

FROM us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:1.0.1
ENV PYTHONDONTWRITEBYTECODE=1

# Rename original so wrapper can import it as _build_vcf_shard_mt_orig
RUN mv /opt/mtSwirl/generate_mtdna_call_mt/Terra/build_vcf_shard_mt.py \
  /opt/mtSwirl/generate_mtdna_call_mt/Terra/_build_vcf_shard_mt_orig.py

# Rename original finalize_mt_with_covdb.py so wrapper can import it as _finalize_mt_with_covdb_orig
RUN mv /opt/mtSwirl/generate_mtdna_call_mt/Terra/finalize_mt_with_covdb.py \
  /opt/mtSwirl/generate_mtdna_call_mt/Terra/_finalize_mt_with_covdb_orig.py

# Rename original add_annotations.py so wrapper can import it as _add_annotations_orig
RUN mv /opt/mtSwirl/generate_mtdna_call_mt/add_annotations.py \
  /opt/mtSwirl/generate_mtdna_call_mt/_add_annotations_orig.py

# Clear any Python bytecode cache that might interfere
RUN find /opt/mtSwirl -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
RUN find /opt/mtSwirl -type f -name "*.pyc" -delete 2>/dev/null || true

# ── VCF pre-localizer wrapper ─────────────────────────────────────────────────
# The wrapper: downloads gs://fake-bucket/* VCFs via STORAGE_EMULATOR_HOST,
# patches sys.argv --shard-tsv to a local path, then delegates to the original.
COPY build_vcf_shard_mt_wrapper.py /opt/mtSwirl/generate_mtdna_call_mt/Terra/build_vcf_shard_mt.py
RUN chmod +x /opt/mtSwirl/generate_mtdna_call_mt/Terra/build_vcf_shard_mt.py && \
  ls -lh /opt/mtSwirl/generate_mtdna_call_mt/Terra/build_vcf_shard_mt.py && \
  ls -lh /opt/mtSwirl/generate_mtdna_call_mt/Terra/_build_vcf_shard_mt_orig.py && \
  echo "✓ build_vcf_shard_mt wrapper installed"

# ── BED file pre-localizer wrapper ────────────────────────────────────────────
# Pre-localizes gs://fake-bucket/* BED files before Hail reads them.
COPY finalize_mt_with_covdb_wrapper.py /opt/mtSwirl/generate_mtdna_call_mt/Terra/finalize_mt_with_covdb.py
RUN chmod +x /opt/mtSwirl/generate_mtdna_call_mt/Terra/finalize_mt_with_covdb.py && \
  ls -lh /opt/mtSwirl/generate_mtdna_call_mt/Terra/finalize_mt_with_covdb.py && \
  ls -lh /opt/mtSwirl/generate_mtdna_call_mt/Terra/_finalize_mt_with_covdb_orig.py && \
  echo "✓ finalize_mt_with_covdb wrapper installed"

# ── Resource file pre-localizer wrapper ───────────────────────────────────────
# Pre-localizes gs://fake-bucket/* resource files (phylotree, variant context, etc.)
# before Hail reads them via hl.import_table().
COPY add_annotations_wrapper.py /opt/mtSwirl/generate_mtdna_call_mt/add_annotations.py
RUN chmod +x /opt/mtSwirl/generate_mtdna_call_mt/add_annotations.py && \
  ls -lh /opt/mtSwirl/generate_mtdna_call_mt/add_annotations.py && \
  ls -lh /opt/mtSwirl/generate_mtdna_call_mt/_add_annotations_orig.py && \
  echo "✓ add_annotations wrapper installed"

# ── gcloud stub ───────────────────────────────────────────────────────────────
# Intercepts `storage cp` and `storage objects describe` and redirects to
# fake-gcs-server via curl. All other subcommands fall through to real gcloud.
COPY gcloud_stub.sh /usr/local/bin/gcloud
RUN chmod +x /usr/local/bin/gcloud && \
  ls -lh /usr/local/bin/gcloud && \
  echo "✓ gcloud stub installed"
