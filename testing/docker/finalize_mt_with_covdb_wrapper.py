#!/usr/bin/env python3
"""Local test wrapper for finalize_mt_with_covdb.

This module is installed *in place of* the production finalize_mt_with_covdb.py
inside the local test Docker image (aou-mitochondrial-combine-vcfs-covdb:local-test).
The original is renamed to _finalize_mt_with_covdb_orig.py.

What it does
------------
1. Reads command-line arguments, specifically --artifact-prone-sites-path.
2. If it starts with gs://fake-bucket/, downloads the file from fake-gcs-server
   (using STORAGE_EMULATOR_HOST env var) to a local temp directory using urllib
   (no GCP credentials needed).
3. Replaces the --artifact-prone-sites-path value with the local file path.
4. Delegates to the original script.

This entirely bypasses Hail's Hadoop GCS connector for BED file reads.
"""

import sys
import os
import tempfile
import urllib.parse
import urllib.request
import shutil
import argparse
import logging

_FAKE_BUCKET = "fake-bucket"
_DEFAULT_EMULATOR = "http://host.docker.internal:4443"


def get_emulator_host() -> str:
    """Get emulator host from env, default to host.docker.internal:4443."""
    return os.environ.get("STORAGE_EMULATOR_HOST", _DEFAULT_EMULATOR)


def download_file_from_fake_gcs(gs_path: str, local_dir: str) -> str:
    """Download gs://fake-bucket/* file from fake-gcs-server to local_dir."""
    # Parse gs://bucket/object path
    if not gs_path.startswith("gs://"):
        raise ValueError(f"Expected gs:// path, got: {gs_path}")
    
    parts = gs_path[5:].split("/", 1)
    if len(parts) != 2:
        raise ValueError(f"Invalid gs:// path: {gs_path}")
    
    bucket, obj_path = parts
    if bucket != _FAKE_BUCKET:
        raise ValueError(f"Only gs://fake-bucket/* is supported for pre-localization")
    
    # Build download URL
    emulator = get_emulator_host()
    obj_path_encoded = urllib.parse.quote(obj_path, safe="")
    url = f"{emulator}/storage/v1/b/{bucket}/o/{obj_path_encoded}?alt=media"
    
    print(f"[finalize_mt_with_covdb_wrapper] Downloading {gs_path} from {url}", file=sys.stderr, flush=True)
    
    # Download to local file
    local_filename = os.path.basename(obj_path)
    local_path = os.path.join(local_dir, local_filename)
    
    try:
        with urllib.request.urlopen(url) as response:
            with open(local_path, "wb") as f:
                shutil.copyfileobj(response, f)
        print(f"[finalize_mt_with_covdb_wrapper] Downloaded to {local_path}", file=sys.stderr, flush=True)
        return local_path
    except Exception as e:
        print(f"[finalize_mt_with_covdb_wrapper] Failed to download {gs_path}: {e}", file=sys.stderr, flush=True)
        raise


def main():
    """Process arguments, download gs:// files locally, then call original script."""
    
    # Create temp directory for downloaded files
    temp_dir = tempfile.mkdtemp(prefix="finalize_mt_covdb_")
    print(f"[finalize_mt_with_covdb_wrapper] Using temp dir: {temp_dir}", file=sys.stderr, flush=True)
    
    # Process sys.argv to replace gs://fake-bucket/* paths with local paths
    new_argv = []
    i = 0
    while i < len(sys.argv):
        arg = sys.argv[i]
        
        if arg == "--artifact-prone-sites-path" and i + 1 < len(sys.argv):
            # Get the path argument
            path_value = sys.argv[i + 1]
            if path_value.startswith("gs://fake-bucket/"):
                # Download the file and replace path
                local_path = download_file_from_fake_gcs(path_value, temp_dir)
                new_argv.append(arg)
                new_argv.append(local_path)
                i += 2
                continue
        
        new_argv.append(arg)
        i += 1
    
    # Replace sys.argv for the original module
    sys.argv = new_argv
    
    print(f"[finalize_mt_with_covdb_wrapper] Modified sys.argv: {sys.argv}", file=sys.stderr, flush=True)
    
    # Execute original as a script
    orig_path = "/opt/mtSwirl/generate_mtdna_call_mt/Terra/_finalize_mt_with_covdb_orig.py"
    
    if not os.path.exists(orig_path):
        print(f"[finalize_mt_with_covdb_wrapper] ERROR: Original module not found at {orig_path}", file=sys.stderr, flush=True)
        sys.exit(1)
    
    # Execute original as a script, providing necessary imports in the namespace
    print(f"[finalize_mt_with_covdb_wrapper] Executing original module: {orig_path}", file=sys.stderr, flush=True)
    try:
        with open(orig_path) as f:
            code = compile(f.read(), orig_path, 'exec')
            # Provide necessary modules in the exec namespace
            exec_globals = {
                '__builtins__': __builtins__,
                '__name__': '__main__',
                '__file__': orig_path,
                'argparse': argparse,
                'logging': logging,
                'os': os,
                'sys': sys,
                'tempfile': tempfile,
                'urllib': urllib,
            }
            exec(code, exec_globals)
    except Exception as e:
        print(f"[finalize_mt_with_covdb_wrapper] FATAL ERROR during execution: {e}", file=sys.stderr, flush=True)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
