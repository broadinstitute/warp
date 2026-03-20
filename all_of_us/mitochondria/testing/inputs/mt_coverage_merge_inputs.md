# mt_coverage_merge_inputs.json Documentation

This file explains the structure and rationale behind the test workflow inputs in `mt_coverage_merge_inputs.json`.

## Input Categories

### 1. Files Localized by miniwdl (Local Paths)

These are main workflow inputs that **miniwdl downloads/localizes before execution**:

```json
"mt_coverage_merge.coverage_tsv": "all_of_us/mitochondria/testing/mocks/coverage.tsv",
"mt_coverage_merge.ancestry_tsv": "all_of_us/mitochondria/testing/mocks/ancestry.tsv",
"mt_coverage_merge.dob_tsv": "all_of_us/mitochondria/testing/mocks/dob.tsv",
"mt_coverage_merge.wgs_median_coverage_tsv": "all_of_us/mitochondria/testing/mocks/wgs_median_coverage.tsv",
"mt_coverage_merge.full_data_tsv": "all_of_us/mitochondria/testing/mocks/full_data.tsv"
```

**Why local paths?**
- Local file paths are faster and avoid requiring GCS credentials
- miniwdl can directly access relative paths without special handling
- This simplifies testing and avoids the cloud-sdk Docker image dependency
- In production on Terra/Cromwell, these would be `gs://` paths and Cromwell would handle localization

### 2. Files Accessed by Hail Docker Image (GCS Paths)

These are paths read **inside the Hail/Spark Docker image**, which expects GCS paths:

```json
"mt_coverage_merge.step3_output_bucket": "gs://fake-bucket",
"mt_coverage_merge.finalize_mt_with_covdb_round1.artifact_prone_sites_path": "gs://fake-bucket/blacklist_sites.hg38.chrM.bed",
"mt_coverage_merge.finalize_mt_with_covdb_round2.artifact_prone_sites_path": "gs://fake-bucket/blacklist_sites.hg38.chrM.bed",
"mt_coverage_merge.finalize_mt_with_covdb_round3.artifact_prone_sites_path": "gs://fake-bucket/blacklist_sites.hg38.chrM.bed"
```

**Why `gs://` paths?**
- Hail uses the Hadoop GCS connector and gcsfs to read GCS paths natively
- These tasks run inside the Hail Docker image, which cannot directly access local file paths
- For testing, the `STORAGE_EMULATOR_HOST=http://host.docker.internal:4443` environment variable redirects all `gs://` requests to fake-gcs-server
- The test Docker image (aou-mitochondrial-combine-vcfs-covdb:local-test) includes pre-localizer wrappers that download these files before Hail reads them
- This approach mirrors production: same input format (gs://), with fake-gcs-server as the backend

## Testing Flow

1. **miniwdl reads inputs from mt_coverage_merge_inputs.json**
2. **miniwdl localizes the main TSV files** (from local paths)
3. **miniwdl runs the workflow**
4. **Hail tasks run inside the test Docker image**
   - Build-time wrappers patch the task commands to redirect gs:// paths
   - Pre-localizers download gs://fake-bucket/* files via STORAGE_EMULATOR_HOST
   - Original Hail code reads the pre-localized files
5. **Results are written to gs://fake-bucket** (intercepted by fake-gcs-server)

## Setting Up Test Data

Test data is uploaded to fake-gcs-server by `all_of_us/mitochondria/testing/setup_fake_gcs.sh`:
- Runs automatically if `run_test.sh` detects missing data
- Uploads coverage files, VCF files, BED files, and mock resource files
- Main TSV files (coverage.tsv, etc.) are NOT uploaded because inputs.json uses local paths for them
