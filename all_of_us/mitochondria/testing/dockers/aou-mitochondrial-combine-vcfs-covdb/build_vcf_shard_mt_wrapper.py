#!/usr/bin/env python3
from __future__ import annotations

# FIRST THING: announce to both stderr AND a debug file
# This tells us if the wrapper is even being executed
import sys
import os

# Determine debug file location - default to current working directory so miniwdl captures it
DEBUG_FILE = os.environ.get("WRAPPER_DEBUG_FILE", "wrapper_startup.log")

def log_debug(msg: str) -> None:
    """Log to stderr and to debug file."""
    print(f"[build_vcf_shard_mt_wrapper] {msg}", file=sys.stderr, flush=True)
    try:
        with open(DEBUG_FILE, "a") as f:
            f.write(f"{msg}\n")
            f.flush()
    except Exception as e:
        print(f"[build_vcf_shard_mt_wrapper] failed to write debug log: {e}", file=sys.stderr, flush=True)

log_debug("STARTED - wrapper script is executing")
log_debug(f"Python {sys.version}")
log_debug(f"sys.argv: {sys.argv}")
log_debug(f"cwd: {os.getcwd()}")

# Verify wrapper environment
ORIG_PATH = "/opt/mtSwirl/generate_mtdna_call_mt/Terra/_build_vcf_shard_mt_orig.py"
if os.path.exists(ORIG_PATH):
    log_debug(f"✓ Original module exists: {ORIG_PATH}")
else:
    log_debug(f"✗ ERROR: Original module NOT found at {ORIG_PATH}")
    log_debug(f"Available files in /opt/mtSwirl/generate_mtdna_call_mt/Terra/:")
    try:
        terra_dir = "/opt/mtSwirl/generate_mtdna_call_mt/Terra"
        if os.path.isdir(terra_dir):
            files = os.listdir(terra_dir)
            for f in files:
                log_debug(f"  - {f}")
        else:
            log_debug(f"  - Terra directory doesn't exist!")
    except Exception as e:
        log_debug(f"  - Failed to list files: {e}")


"""Local test wrapper for build_vcf_shard_mt.

This module is installed *in place of* the production build_vcf_shard_mt.py
inside the local test Docker image (aou-mitochondrial-combine-vcfs-covdb:local-test).
The original is renamed to _build_vcf_shard_mt_orig.py.

What it does
------------
1. Reads the --shard-tsv argument from sys.argv.
2. For every VCF path that starts with ``gs://fake-bucket/``, downloads the
   file from fake-gcs-server (using STORAGE_EMULATOR_HOST env var) to a local
   temp directory using urllib (no GCP credentials needed).
3. Writes a patched shard TSV where those gs:// paths are replaced with the
   local file paths.
4. Replaces the --shard-tsv argument in sys.argv with the patched TSV path.
5. Delegates to the original module's main() / _parse_args().

This entirely bypasses Hail's Hadoop GCS connector for VCF reads, so
STORAGE_EMULATOR_HOST does not need to be wired into Hadoop/Spark config.

For gs:// paths that point to buckets *other* than fake-bucket, the path is
left unchanged (those would fail at Hail level, but that is intentional — the
test only provisions fake-bucket).
"""

import csv
import tempfile
import urllib.parse
import urllib.request

log_debug("dependencies imported successfully")

_FAKE_BUCKET = "fake-bucket"
_DEFAULT_EMULATOR = "http://host.docker.internal:4443"


def _fake_gcs_host() -> str:
    return os.environ.get("STORAGE_EMULATOR_HOST", _DEFAULT_EMULATOR).rstrip("/")


def _localize_vcfs(shard_tsv: str) -> str:
    """Download gs://fake-bucket/* VCFs; return path to a patched shard TSV."""
    host = _fake_gcs_host()
    local_dir = tempfile.mkdtemp(prefix="localized_vcfs_")
    log_debug(f"created temp directory: {local_dir}")

    rows: list[dict] = []
    fieldnames: list[str] = []

    with open(shard_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        log_debug(f"reading shard TSV with fields: {fieldnames}")
        for row_num, row in enumerate(reader, start=1):
            vcf = row.get("vcf", "")
            log_debug(f"row {row_num}: vcf={vcf}")
            gs_prefix = f"gs://{_FAKE_BUCKET}/"
            if vcf.startswith(gs_prefix):
                obj_name = vcf[len(gs_prefix):]
                local_path = os.path.join(local_dir, os.path.basename(obj_name))
                if not os.path.exists(local_path):
                    encoded = urllib.parse.quote(obj_name, safe="")
                    # Use correct fake-gcs-server download endpoint
                    url = f"{host}/storage/v1/b/{_FAKE_BUCKET}/o/{encoded}?alt=media"
                    log_debug(f"downloading {vcf} → {local_path}")
                    log_debug(f"download URL: {url}")
                    try:
                        urllib.request.urlretrieve(url, local_path)
                        log_debug("download OK")
                    except Exception as exc:
                        log_debug(f"download FAILED: {exc}")
                        raise
                else:
                    log_debug(f"file already exists: {local_path}")
                row["vcf"] = local_path
                log_debug(f"row {row_num}: patched vcf={row['vcf']}")
            rows.append(row)

    patched = os.path.join(local_dir, "patched_shard.tsv")
    with open(patched, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    log_debug(f"wrote patched shard TSV: {patched}")
    return patched


def _patch_argv() -> None:
    """Replace --shard-tsv value in sys.argv with localized version."""
    log_debug(f"_patch_argv() called with sys.argv: {sys.argv}")
    for i, arg in enumerate(sys.argv):
        if arg == "--shard-tsv" and i + 1 < len(sys.argv):
            original = sys.argv[i + 1]
            log_debug(f"localizing VCFs from shard TSV: {original}")
            try:
                sys.argv[i + 1] = _localize_vcfs(original)
            except Exception as exc:
                log_debug(f"ERROR during VCF localization: {exc}")
                raise
            log_debug(f"shard TSV patched → {sys.argv[i + 1]}")
            return
    # --shard-tsv not found; let original module surface the error
    log_debug("WARNING: --shard-tsv not found in argv, skipping localization")


# Localize before importing the original (which calls _parse_args at import time
# only if __name__ == '__main__', so this is safe).
log_debug(f"about to call _patch_argv(). sys.argv: {sys.argv}")
_patch_argv()

# Ensure /opt/mtSwirl is in Python path for importing the original module.
if "/opt/mtSwirl" not in sys.path:
    sys.path.insert(0, "/opt/mtSwirl")

# Import and run the original implementation.
try:
    log_debug("about to import original module")
    from generate_mtdna_call_mt.Terra._build_vcf_shard_mt_orig import (  # noqa: E402
        _parse_args,
        main,
    )
    log_debug("successfully imported original module")
except ImportError as exc:
    log_debug(f"FATAL: Failed to import original module: {exc}")
    import traceback
    log_debug(f"Traceback: {traceback.format_exc()}")
    raise
except Exception as exc:
    log_debug(f"FATAL: Unexpected error during import: {exc}")
    import traceback
    log_debug(f"Traceback: {traceback.format_exc()}")
    raise

if __name__ == "__main__":
    try:
        log_debug(f"calling main(_parse_args()) with sys.argv: {sys.argv}")
        main(_parse_args())
    except Exception as exc:
        log_debug(f"FATAL: {exc}")
        import traceback
        log_debug(f"Traceback: {traceback.format_exc()}")
        raise
