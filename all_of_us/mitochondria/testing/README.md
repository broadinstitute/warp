# WARP Test Coverage for mt_coverage_merge

This directory contains test files for running `mt_coverage_merge.wdl` locally with miniwdl.

## Prerequisites

### fake-gcs-server
The test workflow references files in Google Cloud Storage (GCS), which are mocked by [fake-gcs-server](https://github.com/fsouza/fake-gcs-server). You must start this before running tests.

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

**That's it!** `run_test.sh` will automatically:
- Check that fake-gcs-server is running (fail with instructions if not)
- Detect if test data is missing from fake-gcs-server and run `all_of_us/mitochondria/testing/setup_fake_gcs.sh` if needed
- Build the test Docker image from `all_of_us/mitochondria/testing/dockers/aou-mitochondrial-combine-vcfs-covdb/`
- Patch the WDL to use the local test image

## Files

### Directory Structure

```
testing/
├── mocks/                           # Mock data files for testing
│   ├── coverage.tsv                # Metadata coverage table
│   ├── coverage_s001.tsv          # Per-base mtDNA coverage (sample s001)
│   ├── coverage_s002.tsv          # Per-base mtDNA coverage (sample s002)
│   ├── ancestry.tsv               # Ancestry predictions
│   ├── dob.tsv                    # Date of birth data
│   ├── wgs_median_coverage.tsv    # WGS coverage
│   ├── full_data.tsv              # Sample metadata and VCF references
│   ├── s001.vcf                   # Minimal mtDNA VCF (sample s001)
│   ├── s002.vcf                   # Minimal mtDNA VCF (sample s002)
│   ├── blacklist_sites.hg38.chrM.bed # Artifact-prone regions BED
│   ├── phylotree_mock.tsv         # Mock phylogenetic tree data
│   ├── variant_context_mock.tsv   # Mock variant context data
│   ├── pon_mt_trna_mock.tsv       # Mock tRNA PON data
│   └── mitotip_mock.tsv           # Mock MitoTIP pathogenicity scores
├── dockers/                        # Docker test images
│   └── aou-mitochondrial-combine-vcfs-covdb/ # Docker test image
├── inputs/                         # Input data for tests
│   ├── mt_coverage_merge_inputs.json  # WDL input configuration
│   └── mt_coverage_merge_inputs.md   # Input parameter documentation
├── setup_fake_gcs.sh              # Script to setup fake-gcs test data
└── README.md                       # This file
```

### Input Data (TSVs) in `mocks/` subdirectory
- **full_data.tsv** - Main metadata table with sample IDs and references to other data
  - Located in `mocks/`
  - Columns: research_id, contamination, coverage_metrics, final_base_level_coverage_metrics, major_haplogroup, mean_coverage, median_coverage, mtdna_consensus_overlaps, final_vcf
  - `coverage_metrics`: GCS path to per-base mtDNA coverage TSV for the sample
  - `final_base_level_coverage_metrics`: GCS path to per-base mtDNA coverage TSV (renamed to `coverage` in the workflow, which is the default column read by `build_coverage_db`)
  - Both coverage columns point to the same file in this test setup

- **coverage.tsv** - Coverage data keyed by research_id (in `mocks/`)
  - Columns: research_id, mean_coverage, biosample_collection_date, verify_bam_id2_contamination

- **ancestry.tsv** - Genetic ancestry predictions (in `mocks/`)
  - Columns: research_id, ancestry_pred

- **dob.tsv** - Date of birth data (in `mocks/`)
  - Columns: research_id, date_of_birth

- **wgs_median_coverage.tsv** - Whole-genome sequencing coverage (in `mocks/`)
  - Columns: research_id, median_coverage

### Position-Level Coverage Data (in `mocks/` subdirectory)
- **coverage_s001.tsv** - mtDNA position-level coverage for sample s001
  - Columns: position, coverage
  - Referenced in full_data.tsv as `gs://fake-bucket/coverage_s001.tsv`
  
- **coverage_s002.tsv** - mtDNA position-level coverage for sample s002
  - Columns: position, coverage
  - Referenced in full_data.tsv as `gs://fake-bucket/coverage_s002.tsv`

### Mock VCF Files (in `mocks/` subdirectory)
- **s001.vcf** - Minimal GRCh38 (`chrM`) VCF for sample s001
  - Includes all FORMAT fields required by `build_vcf_shard_mt` with `include_extra_v2_fields=true`:
    `GT`, `AD`, `DP`, `AF`, `MQ`, `TLOD`, `FT`, `F1R2`, `F2R1`, `MPOS`
  - 2 variants on chrM (pos 750, 1438)
  - Referenced in full_data.tsv as `gs://fake-bucket/s001.vcf`

- **s002.vcf** - Minimal GRCh38 (`chrM`) VCF for sample s002
  - Same FORMAT fields as s001.vcf; 2 variants on chrM (pos 750, 4769)
  - Referenced in full_data.tsv as `gs://fake-bucket/s002.vcf`

**Note:** These local files must be uploaded to fake-gcs-server before running tests.

### Configuration & Resources
- **mt_coverage_merge_inputs.json** - WDL workflow inputs with two categories:

  **Local paths (miniwdl-localized):**
  - `coverage_tsv`, `ancestry_tsv`, `dob_tsv`, `wgs_median_coverage_tsv`, `full_data_tsv`
  - These are workflow metadata inputs that miniwdl downloads/localizes before execution
  - Using local file paths avoids requiring the cloud-sdk Docker image
  - In production on Terra, these would be gs:// paths (and Cromwell would handle localization)
  
  **GCS paths (Hail-accessed):**
  - `step3_output_bucket` - Output target for Hail tasks
  - `artifact_prone_sites_path` - BED files accessed by Hail/Spark tasks
  - These paths are read inside the Hail Docker image, which expects gs:// URIs
  - For testing, the `STORAGE_EMULATOR_HOST` env var redirects gs:// requests to fake-gcs-server
  - The test Docker image pre-localizers download these files for Hail

- **blacklist_sites.hg38.chrM.bed** - BED file with blacklisted genomic regions (uploaded to fake-gcs, accessed via gs:// by Hail)
- **phylotree_mock.tsv** - Mock phylogenetic tree data (for annotation, not yet used in tests)
- **pon_mt_trna_mock.tsv** - Panel of normals tRNA variants (for annotation, not yet used in tests)
- **variant_context_mock.tsv** - Mock variant context data (for annotation, not yet used in tests)
- **mitotip_mock.tsv** - Mock MitoTip pathogenicity scores (not yet used in tests)

## How It Works

### Setup Flow

1. **Start fake-gcs-server** (see Prerequisites above)
2. **Run `run_test.sh`** from the repository root on your host machine

The script then:
1. **Verifies fake-gcs-server** is running at http://host.docker.internal:4443
2. **Checks for test data** by querying fake-gcs-server for required files
3. **Runs automatic setup** if test data is missing (runs `testing/setup_fake_gcs.sh`)
4. **Builds the test Docker image** from `testing/aou-mitochondrial-combine-vcfs-covdb/`
5. **Patches the WDL** to use local test images (only for images with test versions)
6. **Runs miniwdl** with the patched workflow

## Running Tests

After starting fake-gcs-server, run from the repository root on your host machine:

```bash
# Run with default test inputs (mt_coverage_merge.wdl + all_of_us/mitochondria/testing/inputs/mt_coverage_merge_inputs.json)
./run_test.sh

# Run with custom workflow and test inputs
./run_test.sh . all_of_us/mitochondria/your_pipeline.wdl all_of_us/mitochondria/testing/your_custom_inputs.json

# Skip call cache cleanup (calls cache is always enabled)
SKIP_CLEANUP=true ./run_test.sh

# Force fake-gcs test data reset
RESET_FAKE_GCS=true ./run_test.sh

# Skip Docker image build (use existing test image)
SKIP_BUILD=true ./run_test.sh
```

## Test Coverage

Verified working tasks:
- ✅ `process_tsv_files` - TSV merging and data transformation
- ✅ `make_vcf_shards_from_tsv` - VCF shard creation from full_data.tsv
- ✅ `annotate_coverage` - Coverage TSV annotation via gcsfs (with fake-gcs support)
- ✅ `build_vcf_shard_mt` - Hail VCF shard building (via local test Docker image with pre-localizer)

Currently failing / not yet reached:
- ⏳ `merge_mt_shards` - MatrixTable merge (requires real Hail/Spark setup)
- ⏳ `finalize_mt_with_covdb` - Coverage database finalization
- ⏳ `combine_vcfs_and_homref_from_covdb` - VCF combining and homref addition
- ⏳ `add_annotations` - Resource file annotation (phylotree, variant context, etc.)

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

`run_test.sh` now automatically builds this test image from `all_of_us/mitochondria/testing/dockers/aou-mitochondrial-combine-vcfs-covdb/` and
shadow-tags it as the production image before running miniwdl. Only Docker images with test versions
available (indicated by a corresponding directory under `all_of_us/mitochondria/testing/dockers/`) are patched during WDL processing.
Two hooks are applied:

1. **VCF pre-localizer** (`all_of_us/mitochondria/testing/dockers/aou-mitochondrial-combine-vcfs-covdb/build_vcf_shard_mt_wrapper.py`)  
   Replaces `build_vcf_shard_mt.py` in the image. Before Hail runs, downloads every
   `gs://fake-bucket/*` VCF from fake-gcs-server to local disk (via `urllib` using
   `STORAGE_EMULATOR_HOST`), writes a patched shard TSV with local paths, then delegates
   to the original `build_vcf_shard_mt.py`. Hail reads local files — no GCS connector needed.

2. **`gcloud` stub** (`all_of_us/mitochondria/testing/dockers/aou-mitochondrial-combine-vcfs-covdb/gcloud_stub.sh`)  
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
