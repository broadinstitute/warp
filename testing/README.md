# WARP Test Coverage for mt_coverage_merge

This directory contains test files for running `mt_coverage_merge.wdl` locally with miniwdl.

## Prerequisites

### fake-gcs-server
The test workflow references files in Google Cloud Storage (GCS), which are mocked by [fake-gcs-server](https://github.com/fsouza/fake-gcs-server).

**Start fake-gcs-server before running tests:**

```bash
docker run -d -p 4443:4443 --name fake-gcs \
  fsouza/fake-gcs-server -scheme http \
  -external-url http://host.docker.internal:4443 \
  -public-host host.docker.internal:4443
```

> **Important:** Use `host.docker.internal` (not `localhost`) for `-external-url` and `-public-host`.
> Docker Swarm task containers resolve `localhost` as themselves, not the Mac host, so
> gcsfs would fail to reach fake-gcs-server on any redirect. `host.docker.internal` is
> Docker Desktop's DNS alias for the host machine and works from both host and containers.

The test script will automatically check for fake-gcs-server availability at startup and fail with instructions if it's not running.

**Setup test data in fake-gcs-server:**

```bash
testing/setup_fake_gcs.sh ~/codes/misc/warp
```

## Files

### Input Data (TSVs)
- **full_data.tsv** - Main metadata table with sample IDs and references to other data
  - Columns: research_id, contamination, coverage_metrics, final_base_level_coverage_metrics, major_haplogroup, mean_coverage, median_coverage, mtdna_consensus_overlaps, final_vcf
  - `coverage_metrics`: GCS path to per-base mtDNA coverage TSV for the sample
  - `final_base_level_coverage_metrics`: GCS path to per-base mtDNA coverage TSV (renamed to `coverage` in the workflow, which is the default column read by `build_coverage_db`)
  - Both coverage columns point to the same file in this test setup

- **coverage.tsv** - Coverage data keyed by research_id
  - Columns: research_id, mean_coverage, biosample_collection_date, verify_bam_id2_contamination

- **ancestry.tsv** - Genetic ancestry predictions
  - Columns: research_id, ancestry_pred

- **dob.tsv** - Date of birth data
  - Columns: research_id, date_of_birth

- **wgs_median_coverage.tsv** - Whole-genome sequencing coverage
  - Columns: research_id, median_coverage

### Position-Level Coverage Data
- **coverage_s001.tsv** - mtDNA position-level coverage for sample s001
  - Columns: position, coverage
  - Referenced in full_data.tsv as `gs://fake-bucket/coverage_s001.tsv`
  
- **coverage_s002.tsv** - mtDNA position-level coverage for sample s002
  - Columns: position, coverage
  - Referenced in full_data.tsv as `gs://fake-bucket/coverage_s002.tsv`

### Mock VCF Files
- **s001.vcf** - Minimal GRCh38 (`chrM`) VCF for sample s001
  - Includes all FORMAT fields required by `build_vcf_shard_mt` with `include_extra_v2_fields=true`:
    `GT`, `AD`, `DP`, `AF`, `MQ`, `TLOD`, `FT`, `F1R2`, `F2R1`, `MPOS`
  - 2 variants on chrM (pos 750, 1438)
  - Referenced in full_data.tsv as `gs://fake-bucket/s001.vcf`

- **s002.vcf** - Minimal GRCh38 (`chrM`) VCF for sample s002
  - Same FORMAT fields as s001.vcf; 2 variants on chrM (pos 750, 4769)
  - Referenced in full_data.tsv as `gs://fake-bucket/s002.vcf`

**Note:** These local files must be uploaded to fake-gcs-server before running tests.

### Configuration
- **inputs.json** - WDL workflow inputs pointing to test TSV files
  - Uses GCS paths (e.g., `gs://fake-bucket/...`) that are resolved to fake-gcs-server

## Setting Up fake-gcs-server

### 1. Start the server

```bash
docker run -d -p 4443:4443 --name fake-gcs \
  fsouza/fake-gcs-server -scheme http \
  -external-url http://host.docker.internal:4443 \
  -public-host host.docker.internal:4443
```

### 2. Setup test data in fake-gcs-server

Run the setup script to create the test bucket and upload coverage files:

```bash
testing/setup_fake_gcs.sh ~/codes/misc/warp
```

This script:
- Verifies fake-gcs-server is running
- Creates the `fake-bucket` bucket
- Uploads `coverage_s001.tsv` and `coverage_s002.tsv` to the bucket

After setup completes, you can run the workflow tests.

## Running Tests

From the repository root on your host machine:

```bash
# Run with default test inputs
~/codes/misc/warp/run_test.sh ~/codes/misc/warp

# Run with specific test inputs
~/codes/misc/warp/run_test.sh ~/codes/misc/warp all_of_us/mitochondria/mt_coverage_merge.wdl testing/inputs.json

# Skip cache cleanup
SKIP_CLEANUP=true ~/codes/misc/warp/run_test.sh ~/codes/misc/warp

# Get help
~/codes/misc/warp/run_test.sh -h
```

## Test Coverage

Verified working tasks:
- ✅ `process_tsv_files` - TSV merging and data transformation
- ✅ `make_vcf_shards_from_tsv` - VCF shard creation
- ✅ `annotate_coverage` - Fixed: use `host.docker.internal` for fake-gcs `-public-host`

Currently failing / not yet reached:
- ⏳ `build_vcf_shard_mt` - local test Docker image pre-localizes VCFs + stubs gcloud (in progress)
- ⏳ `merge_mt_shards`, `finalize_mt_with_covdb` - blocked on `build_vcf_shard_mt`
- ⏳ `combine_vcfs_and_homref_from_covdb`, `add_annotations` - downstream of above

Not yet reached (downstream):
- ⏳ `build_vcf_shard_mt`, `merge_mt_shards`, `finalize_mt_with_covdb`
- ⏳ `combine_vcfs_and_homref_from_covdb`, `add_annotations`

## Known Issues

### `annotate_coverage` — gcsfs 401 Unauthorized

**Error:**
```
gcsfs.retry.HttpError: Anonymous caller does not have storage.objects.get access ...
Permission 'storage.objects.get' denied on resource (or it may not exist)., 401
```

**Root cause:**  
Two issues combined:

1. **`-public-host localhost:4443`** — fake-gcs-server was started with `localhost` as the public
   hostname. When gcsfs (inside a Docker Swarm service container) called `_info()` on a
   `gs://` path, gcsfs issued a list-objects request. fake-gcs returned responses containing
   `localhost:4443` in URLs, which the container resolved as itself — not the Mac host — causing
   `FileNotFoundError` or a 401 from the real GCS.

2. **`STORAGE_EMULATOR_HOST` not previously reaching containers** — an earlier run had a
   different fake-gcs startup (scheme differences), and the env var path was unverified.

**Fix applied:**  
Restart fake-gcs with `host.docker.internal` as both `-external-url` and `-public-host`:
```bash
docker rm -f fake-gcs
docker run -d -p 4443:4443 --name fake-gcs \
  fsouza/fake-gcs-server -scheme http \
  -external-url http://host.docker.internal:4443 \
  -public-host host.docker.internal:4443
```
`run_test.sh` already passes `--env STORAGE_EMULATOR_HOST=http://host.docker.internal:4443`
which injects the variable into every Swarm task container. With the public-host fix,
gcsfs routes all requests (including list-objects) to fake-gcs-server correctly.

### `build_vcf_shard_mt` — Hail Hadoop GCS connector

`build_vcf_shard_mt` uses **Hail** (via `gcs-connector-hadoop2-2.2.7.jar`) to read
`gs://` VCF paths, and uses `gcloud storage cp` to write outputs. Neither the Hadoop
GCS connector nor `gcloud` respect `STORAGE_EMULATOR_HOST`.

**Fix: local test Docker image** (`aou-mitochondrial-combine-vcfs-covdb:local-test`)

`run_test.sh` now automatically builds this test image from `testing/docker/` and
shadow-tags it as the production image before running miniwdl. Two hooks are applied:

1. **VCF pre-localizer** (`testing/docker/build_vcf_shard_mt_wrapper.py`)  
   Replaces `build_vcf_shard_mt.py` in the image. Before Hail runs, downloads every
   `gs://fake-bucket/*` VCF from fake-gcs-server to local disk (via `urllib` using
   `STORAGE_EMULATOR_HOST`), writes a patched shard TSV with local paths, then delegates
   to the original `build_vcf_shard_mt.py`. Hail reads local files — no GCS connector needed.

2. **`gcloud` stub** (`testing/docker/gcloud_stub.sh`)  
   Installed as `/usr/local/bin/gcloud`. Handles:
   - `gcloud storage cp <file> gs://fake-bucket/<obj>` → curl upload to fake-gcs
   - `gcloud storage objects describe gs://fake-bucket/<obj> --format=value(md5Hash)` → curl query

The production image tag is restored to the original ID after each run.

## Notes

The test data is **synthetic** and will hit limitations with tasks that require:
- Real genomics data formats (HDF5 coverage databases, Hail MatrixTables)
- GCS interactions beyond what fake-gcs-server + `STORAGE_EMULATOR_HOST` can emulate
- Hail/Spark MatrixTable operations (`build_vcf_shard_mt` uses `aou-mitochondrial-combine-vcfs-covdb:1.0.1`)

For full end-to-end testing, use Terra/Cromwell with real All of Us data.
