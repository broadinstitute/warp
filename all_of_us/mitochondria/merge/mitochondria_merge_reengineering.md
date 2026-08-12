# mtDNA Coverage Merge Pipeline: Scaling to 535k Samples

## Overview

This document describes the changes made to the `Mitochondria Merge` pipeline to scale it from processing ~1,000 samples to successfully processing **535,000 samples** for the All of Us (AoU) v9 data release. 

The pipeline (`mitochondria_merge.wdl`) takes per-sample mitochondrial DNA (mtDNA) variant calls and coverage data as inputs, merges them into a cohort-wide callset, imputes homoplasmic reference genotypes where coverage supports them, and produces annotated VCF output. Before this rewrite, the pipeline failed — crashing with out-of-memory errors after running for >82 hours — when run on the 50,000 sample subset of the v9 AoU cohort.

---

## The Problem: Why the Original Pipeline Could Not Scale

The original pipeline had three distinct bottlenecks, each of which needed to be addressed independently.

### Problem 1: Coverage annotation — the dense matrix bottleneck

**What the pipeline does:** For each of the ~16,569 mitochondrial genome positions, the pipeline computes per-position summary statistics (mean coverage, median coverage, fraction of samples with >100× and >1000× coverage). These statistics are used to decide whether a sample with a missing genotype call at a given position should be assigned a homoplasmic reference genotype (i.e., the sample's reads covered the position well enough that we can be confident no variant exists).

**What the original code did:** Used a [Hail](https://hail.is) MatrixTable (a distributed data structure built on Apache Spark) to represent coverage as a dense **(positions × samples)** matrix. The matrix was read from per-sample coverage TSV files (one per sample), merged into a single matrix, and then row aggregations were run to produce per-position statistics.

**Why it failed at 535k samples:** Hail MatrixTables partition primarily by *rows* (positions). For mtDNA there are only ≤16,569 rows, which means Spark had very few partitions to work with — parallelism was fundamentally limited. More importantly, the full coverage matrix at 535k samples occupied enormous memory on the Spark driver. Operations such as computing exact median coverage per position (a column-wise operation over 535k values) forced the driver to coordinate massive amounts of data. The result was an out-of-memory crash after ~82 hours of runtime.

**Why this was structurally hard to fix within Hail:** In a Hail MatrixTable, columns are always samples and rows are always variants — this is the fixed convention. Hail's distributed execution model partitions data by *rows*, so it is designed and optimized for datasets with many variants (rows) and a comparatively smaller number of samples (columns). For a typical whole-genome cohort the ratio is favorable: millions of variants × thousands of samples. The mtDNA coverage matrix inverts this ratio: only 16,569 position-rows but 535k sample-columns. With so few rows, there are too few partitions for Spark to exploit parallelism, and column-wise operations (such as computing per-position median across 535k samples) require the driver to coordinate enormous amounts of data with no way to distribute the work further. No tuning of Spark parameters could fix this fundamental mismatch between Hail's row-partitioned execution model and the shape of the coverage matrix.

### Problem 2: VCF ingestion — the single-driver bottleneck

**What the pipeline does:** Reads ~535,000 per-sample VCF files from Google Cloud Storage and merges them into a single combined Hail MatrixTable (the "combined VCF MT"), which is the foundation of the annotated cohort callset.

**What the original code did:** Attempted to read all 535k VCF files into a single Hail `import_vcf` call, using a hierarchical merge strategy (`multi_way_union_mts`). The driver process was responsible for coordinating all file opens and orchestrating the multi-way merge in a single monolithic task.

**Why it failed at 535k samples:** Even though each VCF file is small (mtDNA only), at 535k files the per-file overhead — GCS file open latency, metadata coordination, Spark driver memory pressure — became overwhelming. The driver had to track hundreds of thousands of files simultaneously, which exhausted driver memory and produced OOM failures. The problem was not the data volume per se, but the coordination overhead of operating on all 535k files in a single task.

### Problem 3: The finalize step — adapting hom-ref imputation to a new coverage format

**What the pipeline does:** After the combined MT is built, a "finalize" step imputes homoplasmic-reference genotypes: for each sample × position entry where the genotype is missing, if the sample's coverage at that position exceeds a threshold, the genotype is set to hom-ref.

**What the original code did:** The original `determine_hom_refs` function was a single, clean Hail operation. It joined the coverage MatrixTable (produced by `annotate_coverage`) to the combined VCF MT using a standard Hail MT-to-MT join (`coverages[mt.locus, mt.s].coverage`), then applied the hom-ref imputation logic in one `annotate_entries` pass — a single scan over the entire entry space. This was correct and efficient for the original scale.

**Why it needed to change:** Solution 1 replaced the Hail coverage MatrixTable with an HDF5 file (`coverage.h5`). HDF5 is not a Hail object and cannot be joined with a Hail MT-to-MT join expression. Coverage values must instead be read on the Python driver using numpy, then broadcast to Spark workers. Because the full coverage matrix at 535k samples is too large to broadcast all at once, the new implementation must process positions in blocks — reading one block of positions at a time from HDF5 and annotating only those rows.

**The transitional bug introduced during the refactor:** The first implementation of `determine_hom_refs_from_covdb` correctly read HDF5 in blocks, but called `mt.annotate_entries(...)` on the **full, unfiltered MT** for each block — broadcasting a block's coverage slice and updating every entry, with a conditional to no-op entries outside the block. This reproduced the logical result of the original single-pass but scanned all entries ~16–20 times (once per block), multiplying the work by the number of position blocks. With 535k samples and a matrix that is ~99.93% missing, the total entry space is still enormous and scanning it repeatedly was impractical.

**Why this was a regression from the original:** The original single-pass `annotate_entries` was efficient precisely because the coverage MT join was expressed as a Hail lazy expression — Spark evaluated it in one distributed scan with no repeated work. The block-loop approach, applied naively to the full MT, traded that single-pass efficiency for the flexibility of HDF5 I/O, but introduced the repeated-scan cost in doing so.

### Problem 4: Hail writing to `/tmp` (RAM-backed tmpfs)

**What happened:** During `add_annotations`, Hail's `export_vcf` call in CONCATENATED mode writes VCF partition shards to a temporary directory before concatenating them. The `hl.init()` call was placed at module scope in `add_annotations.py` without specifying `tmp_dir`, so Hail defaulted to `/tmp` — which on GCP VMs is a RAM-backed filesystem (tmpfs), limited to ~50% of available RAM. With 535k samples generating many VCF partition shards, the RAM-backed `/tmp` was exhausted, producing `IOException: No space left on device` despite 2 TB of disk space being free.

---

## The Solution: A Redesigned Pipeline Architecture

The v2 pipeline addresses each bottleneck with a targeted, independently testable fix. The guiding principles were:
- **Preserve exact outputs**: all per-position statistics, hom-ref assignments, and cohort-wide variant metrics must match v1 exactly
- **Use the simplest tool that works**: avoid Hail/Spark where a simpler approach suffices
- **Make the pipeline retryable**: WDL tasks are isolated so failed steps can be rerun without restarting everything

### Solution 1: Replace the Hail coverage MT with an HDF5 coverage database

**New component:** `build_coverage_db` task / `build_coverage_db.py` script

**What changed:** The `annotate_coverage` task was replaced with `build_coverage_db`, which streams per-sample coverage TSVs directly into an on-disk **HDF5** (`.h5`) file — entirely bypassing Hail and Spark.

**Why HDF5:** The pipeline runs on a single GCP VM with local-attached disk (Cromwell's local disk). HDF5 is a good fit because it:
- Stores the entire coverage matrix as a single file artifact (simpler WDL I/O and tarball handling)
- Supports partial reads: reading a block of positions across all samples is a single array slice (`coverage[sample_indices, pos_block]`)
- Supports chunking and compression for efficient storage
- Requires only `numpy` and `h5py` — no Spark cluster

**How it works:**
1. The script reads the processed sample TSV to get (sample_id, coverage_TSV_path) pairs
2. It opens each per-sample coverage TSV and writes coverage values into an HDF5 dataset `/coverage` with shape `(n_samples, n_positions)`, dtype `uint16`
3. Per-position statistics are computed exactly:
   - **Mean**: exact floating-point sum divided by N
   - **Median**: computed blockwise using `numpy.partition` (for even N, averages the two middle values and casts to `int32` to match the Hail v1 schema)
   - **over_100 / over_1000**: exact counts, no approximation

**Memory profile:** Processing a block of 256 positions across 535k samples requires only `535,000 × 256 × 2 bytes ≈ 274 MB` of RAM. This is small enough to fit comfortably on any of the pipeline's VMs.

**Output:** `coverage.h5` (plus `coverage_summary.tsv`), packaged as a tarball (`output_db.tar.gz`) passed to downstream tasks.

---

### Solution 2: Sharded VCF ingestion with a merge tree

**New components:** `make_vcf_shards_from_tsv`, `build_vcf_shard_mt`, `make_mt_merge_groups`, `merge_mt_shards` tasks

**What changed:** Instead of one monolithic "import all 535k VCFs" task, the pipeline now uses a **scatter–reduce tree**:

1. **Shard**: `make_vcf_shards_from_tsv` divides the sample list into groups of 2,500 samples each, producing ~214 shard TSV files
2. **Scatter**: `build_vcf_shard_mt` is called in parallel on each shard, importing 2,500 VCFs and writing a shard MatrixTable tarball. Each task runs independently and is individually retryable
3. **Reduce (merge tree)**: `make_mt_merge_groups` + `merge_mt_shards` merge shard MTs in multiple rounds with a fan-in of 10:
   - Round 1: 214 shards → ~22 merged MTs
   - Round 2: 22 merged MTs → ~3 merged MTs  
   - Round 3 (if needed): 3 merged MTs → 1 final merged MT

**Why this scales:**
- Each `build_vcf_shard_mt` task sees only 2,500 VCF files — well within Spark driver capacity
- All shard build tasks run in parallel in Cromwell, limited only by Terra's scatter concurrency
- The merge tree does hierarchical `multi_way_union_mts` (the same merge logic as v1) on manageable input counts at each level
- Intermediate results are persisted as MT tarballs in GCS, enabling call-caching and restart from any point

**WDL artifact handoff:** Intermediate MatrixTables are exchanged between tasks as `.tar.gz` archives. Each task untars to local disk, operates locally, and packs the result for the next task. This avoids ephemeral local path dependencies across WDL task boundaries.

---

### Solution 3: Filter-then-annotate per block, with sample sharding

**New components:** `shard_mt_by_samples`, `finalize_mt_with_covdb`, `union_mt_shards` tasks

**The fix for the transitional bug:** The repeated full-matrix scan was eliminated by adding `mt.filter_rows(mt.__block == block_id)` *before* calling `annotate_entries` for each block. This restricts each annotation pass to only the rows (positions) in the current block. The resulting per-block MTs are then unioned back together with `union_rows_tree`, a fan-in merge to keep the Spark DAG size manageable. The total work is now proportional to the data size — each entry is touched exactly once across all blocks.

**The additional scaling fix — sample sharding:** Even with the per-block filter in place, running `determine_hom_refs_from_covdb` on the full 535k-sample MT in a single task was still too large for one VM. The finalize step was therefore also sharded by samples:

1. **`shard_mt_by_samples`**: splits the combined MT into shards of ~25,000 samples each (~22 shards for 535k samples)
2. **`finalize_mt_with_covdb`** (scattered): each shard runs independently — applying `determine_hom_refs_from_covdb` (with `filter_rows` per block), `remove_genotype_filters`, and `apply_mito_artifact_filter`
3. **`union_mt_shards`**: combines the finalized sample shards back into a single cohort-wide MT via `union_cols`

**The finalize algorithm (covdb-based hom-ref) in the final form:** For each position block within a sample shard:
- `mt.filter_rows(mt.__block == block_id)` — isolate only the rows for this block
- Read `coverage[shard_sample_indices, pos_block_indices]` from HDF5 — only the samples in this shard, only the positions in this block
- Broadcast the coverage slice and annotate entries: if `HL` is missing and `DP > threshold`, set `HL=0.0`, `FT={"PASS"}`, `DP=coverage`
- Checkpoint the block result
- Union all per-block MTs at the end

**Why this scales:**
- Each entry is touched exactly once (filter-then-annotate, not annotate-with-conditional)
- Each finalize shard operates on ~25,000 samples × 16,569 positions — dramatically smaller than 535k × 16,569
- The HDF5 read is a targeted slice `coverage[shard_indices, pos_block_indices]` — not a full file scan
- All shard tasks run in parallel in Cromwell

---

### Solution 4: Minor fixes to `add_annotations.py`

The `add_annotations` step computes cohort-wide variant statistics and exports the annotated VCF. The plan considered a substantial map-reduce redesign of this step, but the existing logic was found to be scalable enough for 535k samples as-is. Only two targeted fixes were needed.

#### Fix 4a: Stop Hail from writing to RAM-backed `/tmp`

**What changed:** `hl.init()` was moved from module scope into the `main()` function, called only after `args.temp_dir` is parsed from the command line. The correct temp directory (on the Cromwell data disk, passed via `--temp-dir ./tmp` in the WDL task command) is now explicitly set for all three of Hail's temp-dir parameters:

```python
hl.init(
    log='annotations_logging.log',
    tmp_dir=f"file://{os.path.abspath(temp_dir)}",
    local_tmpdir=f"file://{os.path.abspath(temp_dir)}",
    spark_conf={"spark.local.dir": os.path.abspath(temp_dir)},
)
```

Additionally, `hl._set_flags(no_whole_stage_codegen="1")` — which was at module scope and triggered Hail's lazy auto-initialization (with default `/tmp`) before `main()` could run — was moved to inside `main()`, immediately after `hl.init()`.

**Why this matters:** On a machine with 504 GB RAM, the RAM-backed `/tmp` tmpfs is approximately 250 GB. Hail's `export_vcf` in CONCATENATED mode writes all VCF partition shards to `tmp_dir` before concatenating them into the final output file. With 535k samples producing many large shards, this exhausted the RAM-backed tmpfs and produced `IOException: No space left on device` despite 2 TB of free disk space. Directing all Hail temp I/O to the Cromwell data disk (4 TB SSD) entirely avoids this.

#### Fix 4b: Read missing `mt_mean_coverage` values from the coverage DB

**What changed:** The new `fill_missing_mt_mean_coverage_from_covdb` function was added to `add_annotations.py`. This reads per-sample mean coverage directly from the HDF5 coverage database (using the same `covdb_utils` utilities used by the finalize step) for any samples whose `mt_mean_coverage` field is missing from the metadata TSV.

**Why this was needed:** The v1 pipeline computed `mt_mean_coverage` as part of the Hail coverage MatrixTable step, so it was always present in the annotation inputs. In v2, the coverage DB (`coverage.h5`) replaced that Hail MT, but some samples could still be missing this value in the input metadata TSV. Rather than requiring a separate preprocessing step, `add_annotations.py` now fills those missing values on-the-fly from the coverage DB before proceeding with annotation.

---

## Summary of Changes by Pipeline Step

| Step | v1 (original) | v2 (new) | Why changed |
|---|---|---|---|
| **Step 1: TSV preprocessing** | `process_tsv_files` — joins coverage, ancestry, DoB, WGS median coverage TSVs | Unchanged | No scaling issues at this step |
| **Step 2: Coverage stats** | `annotate_coverage` — Hail MatrixTable, Spark aggregation over (positions × samples) | `build_coverage_db` — HDF5 file, numpy blockwise computation | Hail MT is wrong data structure for wide coverage matrix; hits driver OOM |
| **Step 3a: VCF ingestion** | Single monolithic import of all 535k VCFs | `make_vcf_shards_from_tsv` + scattered `build_vcf_shard_mt` (2,500 samples/shard) | Single-driver GCS coordination fails at 535k files |
| **Step 3b: Merging** | Single multi_way_union_mts over all shards | 3-round merge tree with fan-in of 10 | Hierarchical reduction keeps each merge within Spark driver capacity |
| **Step 3c: Finalize** | Single task, single-pass `annotate_entries` joining coverage MT | `shard_mt_by_samples` + scattered `finalize_mt_with_covdb` + `union_mt_shards`; HDF5 read in blocks with `filter_rows` per block | Coverage MT replaced by HDF5 (can't use Hail join); naive block-loop over full MT introduced repeated scans; fixed by filter-then-annotate per block + sample sharding |
| **Step 4: Annotations** | `add_annotations.py` with `hl.init()` at module scope (no `tmp_dir` set) and no coverage DB support | Minor fixes only: `hl.init()` moved to `main()` with explicit temp dirs; `fill_missing_mt_mean_coverage_from_covdb` added; annotation logic otherwise unchanged | Hail defaulted to RAM-backed `/tmp`, causing `No space left on device`; coverage DB format required a small compatibility shim |

---

## New WDL Tasks

| Task | Script | Purpose |
|---|---|---|
| `build_coverage_db` | `coverage_db/build_coverage_db.py` (no Hail) | Build HDF5 coverage database from per-sample coverage TSVs |
| `make_vcf_shards_from_tsv` | Inline Python (pandas) | Partition sample list into shard TSVs of ~2,500 samples |
| `build_vcf_shard_mt` | `Terra/build_vcf_shard_mt.py` | Import one shard of VCFs into a Hail MT tarball |
| `make_mt_merge_groups` | Inline Python | Plan merge-tree groups (fan-in) from a list of MT tarballs |
| `merge_mt_shards` | `Terra/merge_mt_shards.py` | Merge one group of MT tarballs into a single MT tarball |
| `shard_mt_by_samples` | `Terra/shard_mt_by_samples.py` | Split a combined MT into per-sample-shard MT tarballs |
| `finalize_mt_with_covdb` | `Terra/finalize_mt_with_covdb.py` | Apply covdb hom-ref imputation + artifact filter to one shard |
| `union_mt_shards` | `Terra/union_mt_shards.py` | Union finalized sample shards back into a single MT |


---

## New Python Modules

| Module | Location | Purpose |
|---|---|---|
| `build_coverage_db.py` | `coverage_db/` | HDF5 coverage DB builder (no Hail/Spark) |
| `covdb_utils.py` | `mtSwirl_refactor/` | Shared utilities for reading the HDF5 coverage DB; used by both `finalize_mt_with_covdb.py` and `add_annotations.py` |
| `merging_utils.py` | `mtSwirl_refactor/` | Core merging logic: `multi_way_union_mts`, `determine_hom_refs_from_covdb`, `apply_mito_artifact_filter` |
| `merging_constants.py` | `mtSwirl_refactor/` | Constants shared across merging steps (filter names, thresholds, etc.) |

---

## Key Design Decisions and Tradeoffs

### Why HDF5 instead of Parquet/TSV/Hail for coverage?

HDF5 was chosen over Parquet because:
- It supports efficient random-access **column (position) slices** across all rows (samples): `coverage[:, pos_block]` is a single I/O call, not a scan
- It produces a single file artifact, which simplifies WDL I/O (one tarball)
- It works well on POSIX local disk (Cromwell's local attached disk)
- It requires no Spark/Hail infrastructure

Parquet would require reading all rows to extract a column slice. TSV would require loading the entire matrix into RAM. Hail MatrixTables have the row-partition structure mismatch described above.

### Why not use Dataproc (managed Spark cluster) instead?

The pipeline runs on **single-VM Cromwell tasks**, not Dataproc. Dataproc would require a separate cluster provisioning step, significantly more operational complexity, and higher cost for tasks that can be solved without distributed computing. The v2 design gets parallelism from Cromwell's scatter/gather — each WDL task is an independent VM — rather than from within-task Spark parallelism.


### Exact output preservation

All changes were designed to preserve **exact** (bit-for-bit) outputs:
- Coverage statistics: mean, median, over_100, over_1000 computed identically via numpy blockwise reduction
- Hom-ref imputation: identical logic (`HL` missing + `DP > threshold` → set hom-ref), just applied per shard
- Cohort-wide variant metrics (AC/AN/AF, histograms, hap/pop arrays): all computed on the full cohort-wide MT after union, not per shard

---

## What Stayed the Same

- The per-sample single-sample Mutect2 pipeline (`mitochondria_single_sample.wdl`) — no changes
- The VCF merge logic (`multi_way_union_mts` hierarchical merge) — same algorithm, now applied to smaller groups
- The hom-ref imputation semantics — identical logic as v1
- The artifact-prone site filter — same filter, applied once after finalization
- The `add_annotations.py` annotation logic — all cohort-wide statistics, filtering, VCF export, and annotation steps are unchanged; only the Hail initialization location and a small coverage DB compatibility shim were added
- All output schemas — annotated VCF fields, MT structure

---

## Tunable Parameters (v2 Defaults)

| Parameter | Default | Description |
|---|---|---|
| `vcf_merge_shard_size` | 2,500 | Samples per VCF ingest shard |
| `vcf_merge_merge_fanin` | 10 | Merge fan-in per round |
| `vcf_merge_shard_n_partitions` | 192 | Hail partitions for each shard MT |
| `finalize_shard_size` | 25,000 | Samples per finalize shard |
| `finalize_shard_n_partitions` | 256 | Hail partitions per finalize shard MT |
| `finalize_union_n_partitions` | 1,000 | Hail partitions for the final unioned MT |

---

## References

- [Planning document: mt_coverage_merge v2 rewrite plan](mitochondria_merge_v2_plan.md)
- [Planning document: Finalize covdb scaling plan](FINALIZE_COVDB_SCALING_PLAN.md)
- [Pull request #1802: mt merge final version](https://github.com/broadinstitute/warp/pull/1802) — the version used to process AoU v9 data
- [WDL changelog](../mitochondria_merge.changelog.md) — `aou_9.0.0`
