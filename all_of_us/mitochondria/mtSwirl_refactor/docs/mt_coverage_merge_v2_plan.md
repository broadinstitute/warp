# mt_coverage_merge v2 rewrite plan (535k samples, exact outputs)

## Context and goals

The current `mt_coverage_merge.wdl` pipeline fails (OOM after ~82 hours) during `annotate_coverage` when scaling beyond ~50k samples. The root cause is that the pipeline uses Hail MatrixTables to represent *coverage* as a dense (positions × samples) matrix, which is extremely wide for mtDNA (≤16,569 positions but 535k samples). Hail partitions primarily by rows (positions), which limits parallelism and leads to driver-memory pressure when attempting operations like exact median.

**Goals**

* Scale to **535,000 samples** on non-Dataproc single-VM execution.
* Preserve **exact outputs** (no approximate statistics): mean/median coverage per position; thresholds >100 and >1000; homref assignment must match current logic.
* Intermediate artifacts do **not** need to be Hail or VCF.
* Final deliverables must include **annotated VCFs** (filtered and unfiltered), optionally also a Hail MT if still useful.
* Land changes as **small, testable commits**.

---

## Intermediate format decision (single recommendation)

### Recommendation: chunked **HDF5** (via `h5py`) for coverage matrix + Parquet/TSV summaries

We will store the per-sample per-position coverage as an on-disk 2D array in **HDF5**.

**Why HDF5 here**

In this pipeline we’re running on a **single VM** with **local-attached disk** (Cromwell local disk) and we want minimal operational complexity. HDF5 is a good fit because it:

* yields a **single-file artifact** (simplifies WDL I/O and tarball handling)
* supports **chunking + compression** (e.g., `lzf` or `gzip`) for the wide coverage matrix
* supports **partial reads** (read a position block across all samples efficiently)
* has a mature ecosystem and predictable performance on POSIX filesystems

**Coverage DB artifacts**

* `coverage.h5`
  * dataset `/coverage` shape `(n_samples, n_positions)` dtype `uint16` or `uint32`
  * dataset `/sample_id` shape `(n_samples,)` dtype string
  * attribute `position_start=1`, `position_end=16569`
  * attribute `reference="GRCh38"` (positions correspond to MT reference positions as in current coverage file)
* `coverage_summary.tsv`
  * columns: `pos, mean, median, over_100, over_1000`
* Optional (only if needed): `sample_level_coverage.tsv.bgz` (sharded output)

---

## v2 pipeline architecture (single recommendation per step)

### Step 0 — Keep: `subset_data_table`

No functional changes.

**Reason**: It’s not a scaling bottleneck.

**Test**: Compare output row count and sample IDs with the input sample list.



### Step 1 — Keep: `process_tsv_files`

No functional changes planned.

**Reason**: This step is not the scaling bottleneck and the current behavior is sufficient for the v2 rewrite.

**Test**: row count preserved; required columns non-null where expected.

---

### Step 2 — Replace: `annotate_coverage` → `build_coverage_db`

Replace Hail-based `Terra/annotate_coverage.py` with a new script that:

1) Streams per-sample coverage TSVs into `coverage.h5` in batches.
2) Computes per-position metrics exactly:
   * mean (exact)
   * over_100, over_1000 fractions (exact)
   * median (exact) using `numpy.partition` blockwise

**Algorithm outline**

* Build the matrix: write rows as samples, columns as positions.
  * Batch size tuned to RAM (e.g., 2,000–10,000 samples per batch).
* Compute mean / counts in one streaming pass.
* Compute exact median:
  * Process by **position blocks** (e.g., 256 positions at a time).
  * For each block, read `coverage[:, pos_block]` into memory (shape 535k × 256).
  * Compute median per column using `np.partition`.

**Why this scales**

* Memory bounded by one block: `535k * 256 * 2 bytes ≈ 274MB` if `uint16`.
* No construction of massive Hail MT objects.
* Disk I/O is sequential and predictable.

**Outputs**

* `coverage.h5`
* `coverage_summary.tsv`

**Implementation (v2 script)**

* `mtSwirl/generate_mtdna_call_mt/coverage_db/build_coverage_db.py`

**Dependencies**

This script requires Python packages:

* `numpy`
* `h5py`

In local development, you may need to install `h5py` (it is not currently available in this workspace’s Python env).

**Median exactness note**

The v1 coverage MT schema uses `median: int32`. The v2 builder computes an exact median per position. For even $N$, it computes the average of the two middle values and casts to `int32`.
If parity tests indicate Hail rounds differently for even $N$, adjust only that final casting/rounding step.

**Test**

* On 1k samples: run old Hail annotate + new build_coverage_db and compare:
  * per-position mean/median/over_100/over_1000 exactly equal
  * sample-level coverage spot-check equality for a handful of samples/positions

---

### Step 3 — Replace: monolithic combine+impute with **sharded ingest → merge tree → finalize once**

The current Hail `combine_vcfs.py` merges sample VCFs into a combined Hail MatrixTable, then uses per-sample/per-position coverage to impute DP and convert missing calls into homoplasmic reference (HL=0) when coverage is high enough.

In v2, we will keep the *same semantics* for hom-ref / DP imputation, but we will replace the coverage MatrixTable dependency with the HDF5 coverage DB (`coverage.h5`) produced in Step 2.

**Decision (current implementation):** Step 3 produces a **single finalized** Hail MatrixTable (exchanged between tasks as MT tarballs during merging), and runs covdb hom-ref/DP imputation + artifact filtering **once** at the end.

#### Updated Step 3 recommendation for 535k samples: shard ingestion + merge-tree

At ~535k samples, the bottleneck is not VCF *bytes* (VCFs are tiny) but per-file overhead (GCS open latency + driver coordination). The stable approach is:

1) **Shard** the sample list into shard manifests of ~2,500 samples.
2) **Scatter**: import each shard into a multi-sample MT shard.
3) **Reduce**: merge shard MTs in multiple rounds with fan-in 10 (tree reduction).
4) **Finalize once**: after all merges, run covdb hom-ref/DP imputation + artifact-prone site filter exactly once on the final merged MT.

This provides:
* High parallelism and retryability (one shard per task)
* Avoids a single massive “read 535k tiny VCFs” driver bottleneck
* Preserves exact semantics because the final MT is still a single cohort-wide dataset when Step 4 runs

**Recommended defaults** (tunable):
* shard size: **2,500 samples / shard**
* merge fan-in: **10**
* shard MT partitions: **128** (keeps shard MTs reasonably sharded without exploding small files)

**New Terra entrypoints (mtSwirl)**

* `generate_mtdna_call_mt/Terra/make_vcf_shards_from_tsv.py`
  * input: processed TSV
  * output: `vcf_shard_00000.tsv`, … plus `shards.tsv` index
* `generate_mtdna_call_mt/Terra/build_vcf_shard_mt.py`
  * input: one shard manifest TSV
  * output: one shard MT (`vcf_shard_00000.mt`)
* `generate_mtdna_call_mt/Terra/merge_mt_shards.py`
  * input: a TSV of MT paths (fan-in group)
  * output: one merged MT
* `generate_mtdna_call_mt/Terra/finalize_mt_with_covdb.py`
  * input: final merged MT + `coverage.h5`
  * output: finalized MT (covdb homref/DP + artifact filter)

**Artifact handoff (important Terra/WDL detail)**

In the Warp WDL, intermediate Hail MatrixTables are exchanged between tasks as **tarballs** (e.g. `shard_mt.tar.gz`, `merged_mt.tar.gz`). This makes the workflow restartable and ensures Cromwell has concrete `File` artifacts to persist/localize between tasks.

Concretely:

* The shard build step writes a shard MT directory (`*.mt`) and then **tars it**.
* Merge steps consume a list of **MT tarballs** (header `mt_tar`), **untar each to local disk**, and then generate a local TSV (header `mt_path`) that is passed into `merge_mt_shards.py`.
* Finalization consumes the final merged MT **tarball**, untars it locally, then runs covdb hom-ref/DP imputation + artifact-prone site filtering.

This avoids passing ephemeral local filesystem paths between WDL tasks.

**Approach**

* Read the processed TSV to collect per-sample VCF paths.
  * VCFs are **uncompressed** and live in **Google buckets** (e.g. `gs://...`).
  * IMPORTANT: Terra will **not localize** these VCFs because the file paths come from within a TSV and we cannot pass ~535k files as explicit WDL inputs.
  * Therefore, Step 3 must read `gs://` paths directly (similar to how Step 2 reads `gs://` coverage TSVs).

* Merge VCFs into a combined MT using the same hierarchical strategy as v1 (`vcf_merging_and_processing`, `vcf_merging`, `multi_way_union_mts`).

* Replace `determine_hom_refs(mt, coverage_mt_path, minimum_homref_coverage)` with:
  * `determine_hom_refs_from_covdb(mt, coverage_h5_path, minimum_homref_coverage)`
  * Exact parity contract (matches v1):
    1) For entries where `HL` is missing, set `DP = coverage`.
    2) If `HL` is missing and `DP > minimum_homref_coverage`, set `HL=0.0` and `FT={"PASS"}`.
    3) If `HL` is missing and `DP <= minimum_homref_coverage`, set `DP` back to missing.

* Coverage DB lookup contract (from Step 2):
  * `coverage.h5` contains dataset `/coverage` with shape `(n_samples, n_positions)`.
  * `sample_id` is stored in dataset `/sample_id` aligned to the first axis.
  * `pos` is stored in dataset `/pos` aligned to the second axis (1-based positions).
  * For our purposes, sample IDs are in the row index and basepair positions are the columns.

* Output is a single Hail MatrixTable directory (`.mt`) (sharded by partitions on disk) plus any optional lightweight QC exports.

#### Step 3 scaling implementation detail (v2)

To maximize throughput at ~535k samples without requiring a >2TB RAM machine, Step 3 uses a **block-broadcast** covdb annotation strategy:

* Add a stable MT column index (`col_idx`) and map each sample to its covdb row index.
* Sort rows by position and repartition rows so each partition corresponds to a contiguous **position block**.
* For each position block:
  * Read `coverage.h5` for all MT samples but only the columns (positions) in that block.
  * Broadcast *only that block* to executors.
  * Annotate entries with `DP = coverage` for missing HL using `cov_block[offset][col_idx]`.
    * Important: MT rows may be sparse (non-contiguous positions), so the implementation uses a per-block mapping `pos -> offset` rather than assuming contiguous positions.

This avoids materializing a (pos × sample) table and avoids broadcasting a giant `pos -> vector` dictionary. The block size is tunable (e.g., 256–2048 positions) to trade off between broadcast size and number of reads.

**Why this scales**

* We avoid constructing the enormous *coverage* MatrixTable entirely (the original scaling bottleneck).
* Coverage lookup uses block reads from HDF5 (efficient on the single VM local disk).
* The MT output is sharded by Hail partitions, which is compatible with Step 4 and avoids a single monolithic VCF artifact.

**Outputs**

* `combined.mt` (Hail MatrixTable directory; sharded by partitions)

#### Step 3 WDL wiring (current state)

The Warp workflow `warp/all_of_us/mitochondria/mt_coverage_merge.wdl` implements Step 3 as:

* `make_vcf_shards_from_tsv` → shard TSVs (`vcf_shard_*.tsv`)
* `build_vcf_shard_mt` (scatter) → shard MT tarballs
* `make_mt_merge_groups` → merge-group TSVs
* `merge_mt_shards` (scatter, multi-round) → merged MT tarballs
* `finalize_mt_with_covdb` (called once on the deepest merge output) → final MT tarball

#### Step 3 WDL wiring (sharded mode)

`warp/all_of_us/mitochondria/mt_coverage_merge.wdl` now supports a sharded mode controlled by:

* `shard_step3 = true` (default)
* `step3_shard_size = 2500`
* `step3_merge_fanin = 10`

The workflow structure:
* make shard manifests from the processed TSV
* `scatter` build shard MTs
* 3 merge rounds (fan-in 10) to reduce to a single merged MT
* finalize once with covdb + artifact filter and emit the final MT tarball

#### Step 3 Docker image recommendation

Step 3 needs a **Hail/Spark-capable** image *and* the updated mtSwirl code (including the sharded ingest/merge/finalize scripts + `merging_utils.py`).

Recommended: **new dedicated image** (keeps Step 3 changes isolated)

* Dockerfile: `mtSwirl/generate_mtdna_call_mt/Terra/Dockerfile`
* Base image: `us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0`
* Intended tag used by the WDL task:
  * `us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:dev`

Alternative: bake the new script into the existing `aou-mitochondrial-annotate-coverage` image.
This works too, but it’s riskier because that image may be shared by other workflows.

**Logging / observability requirements**

Step 3 must include INFO-level logging suitable for Terra job monitoring.

Minimum recommended logs:

* number of samples discovered in TSV
* progress importing VCFs (e.g., every 1k samples)
* progress through multi-way union stages (already logged in v1 via `multi_way_union_mts`)
* hom-ref annotation phase start/end
* checkpoint/write start/end and output path

**Test**

* On ~1k samples: compare MT entry fields (`DP`, `HL`, `FT`) against v1 output to confirm hom-ref/DP logic matches exactly.

---

### Step 4 — Replace: `add_annotations` with a **two-phase** annotation design

Running `add_annotations.py` independently on sample-sharded callsets will change results, because the script computes **cohort-wide per-variant metrics** using `hl.agg.*` row aggregations.

#### What in `add_annotations.py` is cohort-wide?

Scanning `generate_mtdna_call_mt/add_annotations.py` shows many row aggregations that must see **all samples** to match current output exactly. Examples include (not exhaustive):

* In `generate_expressions()` (row-level metrics):
  * `AC`, `AN`, `AF` (overall)
  * `AC_hom`, `AC_het`, `AF_hom`, `AF_het`
  * `hl_hist` (`hl.agg.hist` of heteroplasmy)
  * `dp_hist_all`, `dp_hist_alt`, `dp_mean`, `mq_mean`, `tlod_mean`
  * `max_hl`
  * haplogroup/population stratified metrics using `hl.agg.group_by(input_mt.hap, ...)` and `hl.agg.group_by(input_mt.pop, ...)`
  * `filtering_allele_frequency(...)` (FAF)
* In `add_quality_histograms()`:
  * global histograms across all variants (`aggregate_rows`) and age distribution across all samples (`aggregate_cols`)
  * per-variant age histograms via `age_hists_expr(..., input_mt.GT, input_mt.age)`
* In `apply_common_low_het_flag()`:
  * `AC_mid_het`, `AF_mid_het`, and downstream `common_low_heteroplasmy`

Because these are computed with `hl.agg` over *all columns*, per-shard execution would compute shard-wide AC/AF/etc instead of cohort-wide values.

#### Schema-driven confirmation: what must be preserved exactly

The schemas in `generate_mtdna_call_mt/schemas.txt` are the target outputs for v2 equivalence:

* **Coverage MT schema** confirms per-position metrics that must match exactly:
  * row: `mean`, `median`, `over_100`, `over_1000` with entry `coverage`
* **Annotated MT schema** confirms the cohort-wide row metrics we must reproduce exactly include:
  * counts: `AN`, `AC_hom`, `AC_het`, `excluded_AC`
  * derived rates: `AF_hom`, `AF_het`
  * summaries: `dp_mean`, `mq_mean`, `tlod_mean`, `max_hl`
  * histograms: `hl_hist`, `dp_hist_all`, `dp_hist_alt`, `age_hist_hom`, `age_hist_het`
  * hap/pop arrays: `hap_*`, `pop_*` plus `hap_hl_hist`, `pop_hl_hist`
  * FAF-related: `hap_faf_hom`, `hapmax_AF_hom`, `hapmax_AF_het`, `faf_hapmax_hom`

This schema review is also strong evidence that for annotated MT cohort-wide row metrics we are **not** relying on per-variant exact medians/quantiles beyond histogram outputs (see “Order-statistics assumption” below).

#### Canonical histogram/binning contract (to prevent Phase 4A/4B drift)

Many cohort-wide metrics in the annotated MT schema are histograms. In v1 they are computed inside Hail using `hl.agg.hist(...)` (directly) and `age_hists_expr(...)` (via gnomAD utilities). To guarantee exact output parity and avoid subtle “off-by-one bin edge” drift when implementing Phase 4A/4B, we define a single canonical binning contract here.

**Rule 1: Phase 4A is responsible for producing the final histogram structs**

Phase 4A must compute (and the reducer must persist) *the exact* `bin_edges`, `bin_freq`, `n_smaller`, and `n_larger` for each per-variant histogram field in the `Annotated_vcf_mt` schema. Phase 4B should not “recompute” histograms; it should only join and possibly reshape them for VCF output.

**Rule 2: Always persist `bin_edges` for any struct-hist field, even if the schema mostly uses bin_freq**

The schema includes histogram structs like `hl_hist`, `dp_hist_all`, `dp_hist_alt`, `age_hist_hom`, `age_hist_het` (each with `bin_edges`, `bin_freq`, `n_smaller`, `n_larger`). Phase 4A output must include these edges explicitly.

**Rule 3: All histogram semantics must match Hail’s `hl.agg.hist`**

Even if we implement Phase 4A without Hail (e.g., in Python/NumPy), the semantics must match Hail’s histogram behavior exactly (including how boundary values are handled and how `n_smaller`/`n_larger` are counted). The lowest-risk approach is to compute histograms with Hail on shards and reduce by summing bin counts.

##### Canonical histogram definitions (as implemented in v1)

Per-variant histograms (row fields):

* `hl_hist`:
  * Code: `generate_expressions()`
  * Definition: `hl.agg.filter(HL > 0, hl.agg.hist(HL, 0, 1, 10))`
  * Output struct fields: `bin_edges`, `bin_freq`, `n_smaller`, `n_larger`

* `dp_hist_all`:
  * Code: `generate_expressions()`
  * Definition: `hl.agg.hist(DP, 0, 2000, 10)`

* `dp_hist_alt`:
  * Code: `generate_expressions()`
  * Definition: `hl.agg.filter(GT.is_non_ref(), hl.agg.hist(DP, 0, 2000, 10))`

* `age_hist_hom`, `age_hist_het`:
  * Code: `add_quality_histograms()`
  * Definition: `age_data = age_hists_expr(True, GT, age)` then `annotate_rows(age_hist_hom=..., age_hist_het=...)`
  * **Important**: `age_hists_expr` is imported from `gnomad.utils.annotations` (not defined in this repo). Treat it as an external contract.
  * v2 guidance:
    * Prefer computing these histograms in Phase 4A with Hail using the same `age_hists_expr` call on shard MTs.
    * Persist the resulting `bin_edges` from the Hail struct; do not hard-code guessed age bins.
    * Reducer merges by elementwise sum of `bin_freq` and sum of `n_smaller`/`n_larger`.

Filter-failure histograms (row fields, arrays of counts only):

* `base_qual_hist`, `map_qual_hist`, `position_hist`, `strand_bias_hist`, `weak_evidence_hist`, `contamination_hist`, `heteroplasmy_below_min_het_threshold_hist`:
  * Code: `add_filter_annotations()` calls `generate_filter_histogram()`
  * Definition: `hl.agg.filter(str(FT).contains(filter_name), hl.agg.hist(HL, 0, 1, 10)).bin_freq`
  * Note: only `bin_freq` is stored (no edges) in the MT schema; nevertheless, the implicit edges are exactly those of `hl.agg.hist(HL, 0, 1, 10)`.

Global histograms (globals in `Annotated_vcf_mt`):

* `dp_hist_all_variants_*`, `mq_hist_all_variants_*`, `tlod_hist_all_variants_*`:
  * Code: `add_quality_histograms()`
  * Definitions:
    * `hist(dp_mean, 0, 4000, 40)`
    * `hist(mq_mean, 0, 80, 40)`
    * `hist(tlod_mean, 0, 40000, 40)`
  * v2 note: these are *aggregate_rows* over variants; because the number of variants is tiny (≤16k), these can be computed in Phase 4B after joining global per-variant means.

* `age_hist_all_samples_*`:
  * Code: `add_quality_histograms()`
  * Definition: `input_mt.aggregate_cols(hl.agg.hist(age, 30, 80, 10))`
  * v2 note: this is a cohort-level histogram across samples, not variants. It is reducible across shards by summing the `bin_freq` and `n_smaller`/`n_larger`.

##### Implementation guardrails

To prevent subtle drift during implementation:

* Phase 4A should output a small JSON/YAML “histogram metadata” blob alongside `global_variant_metrics` that includes each histogram’s `bin_edges` (for struct-hist fields) and the provenance (e.g., `hl.agg.hist(HL,0,1,10)` or `gnomad.utils.annotations.age_hists_expr(True,GT,age)`).
* Phase 4B should validate that bin edge arrays match expectations (cheap assert on first row) before exporting.

#### New recommended approach (exact results): Map/Reduce over shards

We split annotation into **Phase 4A (global metrics)** and **Phase 4B (finalization/export)**.

##### Phase 4A — compute **global cohort-wide variant metrics** once

**Input:** sample-sharded callset shards (VCF or MT shards) + processed metadata TSV.

**Work:**

1) For each sample shard, compute a **partial variant-metrics table** keyed by `(locus, alleles)` containing mergeable components:
   * counts/sums needed for AC/AN/AC_hom/AC_het, dp_mean/mq_mean/tlod_mean (store sums + n)
   * histograms as bin_freq arrays (merge by elementwise sum)
   * group_by hap and pop: store dictionaries or fixed-order arrays of per-group partials
   * any other row-level flags computed with counts (e.g., `AC_mid_het`)
2) Reduce all shard partial tables into a single **global variant metrics table** (≤16k rows).

**Output:** `global_variant_metrics.{ht|parquet|tsv}` keyed by `(locus, alleles)`.

##### Phase 4A sufficient-statistics record spec (source of truth)

Phase 4A must emit a single “sufficient statistics” record per `(locus, alleles)` that lets Phase 4B reconstruct *every* cohort-wide row metric in the `Annotated_vcf_mt` schema deterministically.

We define a reducer-friendly record `VariantMetricsPartial` with the following fields.

**Key**

* `locus: locus<GRCh38>`
* `alleles: array<str>`

**Scalar counts and sums (overall; all samples)**

These are computed on the post-transform entry state used by the current pipeline (after `add_filter_annotations()` updates `FT`/`HL`, and after low-allele-frac + min-het operations).

* `AN: int64`  
  count of defined `HL` (`hl.is_defined(HL)`)
* `AC_hom: int64`  
  count where `HL >= min_hom_threshold` (default 0.95)
* `AC_het: int64`  
  count where `0 < HL < min_hom_threshold`
* `excluded_AC: int64`  
  count where `FT != {'PASS'}` (matches `add_filter_annotations()`)

Means must be reconstructable exactly from integer sums + counts:

* `dp_sum: float64`, `dp_n: int64`  
  sum/count for `DP` (matching `hl.agg.mean(DP)` behavior)
* `mq_sum: float64`, `mq_n: int64`
* `tlod_sum: float64`, `tlod_n: int64`

Max HL:

* `max_hl: float64`  
  max observed `HL` (over all entries, matching `hl.agg.max(HL)`)

**Histograms (overall)**

All histogram bin edges are fixed and must match the code in `add_annotations.py`:

* HL bins: `hl.agg.hist(HL, 0, 1, 10)` (and filtered to `HL > 0` where noted)
* DP bins: `hl.agg.hist(DP, 0, 2000, 10)`
* Age bins (see `age_hists_expr`): must match its hardcoded edges/nbins.

Each histogram is represented by:

* `*_bin_freq: array<int64>`
* `*_n_smaller: int64`
* `*_n_larger: int64`

Required histograms:

* `hl_hist_*`  
  histogram over `HL` with filter `HL > 0`
* `dp_hist_all_*`  
  histogram over `DP` (all entries)
* `dp_hist_alt_*`  
  histogram over `DP` restricted to `GT.is_non_ref()`
* `age_hist_hom_*`  
  histogram over `age` restricted to hom calls (definition must match `age_hists_expr(True, GT, age)`)
* `age_hist_het_*`  
  histogram over `age` restricted to het calls (definition must match `age_hists_expr`)

**Common-low-heteroplasmy sufficient stats**

`common_low_heteroplasmy` is computed in `apply_common_low_het_flag()` from a mid-het AF threshold.
Phase 4A should emit the numerator and denominator so Phase 4B can set the boolean exactly:

* `mid_het_AC: int64`  
  count where `(0 < HL < 0.50) & (FT in {'PASS','low_allele_frac'})`
* `mid_het_AN: int64`  
  count where `is_defined(HL) & (FT in {'PASS','low_allele_frac'})`

Phase 4B reconstructs:

* `common_low_heteroplasmy = (mid_het_AC / mid_het_AN) > 0.001`

**Haplogroup-stratified sufficient stats**

To match `add_annotations_by_hap_and_pop()`, Phase 4B must output arrays aligned to a deterministic `hap_order`. Therefore:

* Phase 4A must emit per-hap partials in a representation that can be reduced without knowing the final global order.
* Recommended: emit as dictionaries keyed by hap string, then Phase 4B standardizes to arrays.

Per hap `h`:

* `hap_AN[h]: int64`
* `hap_AC_hom[h]: int64`
* `hap_AC_het[h]: int64`
* `hap_hl_hist_bin_freq[h]: array<int64>` (HL bin counts)

FAF inputs (computed from counts in Phase 4B):

* `hap_FAF_hom[h]: float64` is computed via `hl.experimental.filtering_allele_frequency(AC_hom, AN, 0.95)`

Phase 4B derived hap metrics:

* `hap_AF_hom[h] = hap_AC_hom[h] / hap_AN[h]`
* `hap_AF_het[h] = hap_AC_het[h] / hap_AN[h]`
* `hap_faf_hom[h] = filtering_allele_frequency(hap_AC_hom[h], hap_AN[h], 0.95)`
* `hapmax_AF_hom`: `hap_order[argmax(hap_AF_hom)]` (must match Hail tie/`unique=True` behavior)
* `hapmax_AF_het`: similarly
* `faf_hapmax_hom = max(hap_faf_hom)`

**Population-stratified sufficient stats**

Same approach as hap (dict keyed by pop string), with Phase 4B producing arrays aligned to `pop_order`:

Per pop `p`:

* `pop_AN[p]: int64`
* `pop_AC_hom[p]: int64`
* `pop_AC_het[p]: int64`
* `pop_hl_hist_bin_freq[p]: array<int64>`

Phase 4B derived:

* `pop_AF_hom[p] = pop_AC_hom[p] / pop_AN[p]`
* `pop_AF_het[p] = pop_AC_het[p] / pop_AN[p]`

**Reducer contract**

The reducer for partial records is strictly associative:

* scalar counts: sum
* sums/n: sum
* max: max
* histogram bin_freq: elementwise sum; `n_smaller/n_larger`: sum
* hap/pop dicts: per-key sum (and elementwise sum for bin arrays)

**Phase 4B join contract**

Phase 4B reads the final `VariantMetricsGlobal` table keyed by `(locus, alleles)` and does:

* Join onto the MT rows
* Emit exactly the row-schema fields in `Annotated_vcf_mt`
* Ensure `hap_order` and `pop_order` in globals match the ordering used to standardize arrays

This table is small (≤16k rows), so Phase 4B can broadcast it efficiently.

**Why this preserves exactness:**

* All cohort-wide metrics are computed as a global reduction, not shard-local approximations.
* All operations are associative/mergeable (sums, counts, histogram bin sums, per-group sums).

##### How Phase 4A and Phase 4B produce a single cohesive output

Phase 4A yields a **small row table** keyed by `(locus, alleles)` (≤16k rows). Phase 4B reads the callset (MT/VCF shards), applies sample-level filters and external resource joins, then **annotates rows by key lookup** from the Phase 4A table before exporting.

Conceptually (Hail MT form):

* `metrics_ht = hl.read_table(global_variant_metrics).key_by('locus','alleles')`
* `mt = mt.annotate_rows(**metrics_ht[mt.row_key])`
* export final outputs (filtered and unfiltered)

This ensures every shard (and any merged output) carries identical cohort-wide INFO/row metrics.

##### Map/Reduce friendliness of cohort-wide metrics (assumption: no exact order-stat metrics)

For the current `add_annotations.py`, the cohort-wide metrics we rely on are built from aggregators that are **mergeable** across shards. That includes:

* **Counts / sums / means**: `count_where`, `mean` (implemented as sum + n internally), and derived rates like AF.
* **Min/max**: `hl.agg.min`, `hl.agg.max`.
* **Histograms**: `hl.agg.hist` (merge by elementwise summing `bin_freq`, and keep shared `bin_edges`).
* **Group-bys** over `hap` and `pop`: merge per-group partials using the same rules above.
* **FAF**: computed from AC/AN (mergeable counts).

**Order-statistics assumption (validated by schema and code paths):**

For cohort-wide *annotated MT* row metrics, we rely only on mergeable aggregators (counts/sums/means/max and histograms). The `Annotated_vcf_mt` schema in `schemas.txt` contains no per-variant median/quantile fields (beyond histogram structs), and the relevant row-metric code paths (`generate_expressions()`, `add_quality_histograms()`, `add_filter_annotations()`, `apply_common_low_het_flag()`, and `add_annotations_by_hap_and_pop()`) do not compute quantiles/medians.

Separately, the **coverage** step does require an exact per-position median (coverage MT row field `median`), which is handled by `build_coverage_db`.

##### Phase 4B — apply external annotations, sample filters, join global metrics, and export

**Input:**

* per-sample shards (same as Phase 4A)
* processed metadata TSV (for `pop`, `age`, etc)
* `global_variant_metrics` (from Phase 4A)

**Work:**

1) Apply sample-level filters (copy number, contamination, overlaps) in a way that matches current script.
2) Apply variant-level resource joins (variant context, phylotree/hap-defining, tRNA predictions, dbSNP/vep if enabled).
3) Join in `global_variant_metrics` by `(locus, alleles)`.
4) Export final deliverables:
   * unfiltered annotated VCF (or annotated hail MT)
   * filtered annotated VCF (or filtered hail MT)

**Output:** final annotated callsets matching current output.

#### Tests

* On 1k samples: run original pipeline and v2 pipeline and compare:
  * global variant metrics fields (AC/AN/AF, histograms, FAF, pop/hap arrays)
  * final annotated VCF content (at least INFO fields; optionally sample FORMAT fields)

---

## WDL changes (v2)

### High-level workflow sketch

1. `subset_data_table` (optional)
2. `process_tsv_files`
3. `build_coverage_db`
4. `make_vcf_shards_from_tsv`
5. `scatter` over shards:
  * `build_vcf_shard_mt`
6. merge rounds (fan-in):
  * `make_mt_merge_groups`
  * `merge_mt_shards`
7. `finalize_mt_with_covdb`

---

## Incremental migration plan (small commits)

### Commit 1 — Add coverage DB builder (no pipeline change)

* Add script: `generate_mtdna_call_mt/coverage_db/build_coverage_db.py`
* Add unit test / smoke test script for 10–100 samples.

**Acceptance**: on 1k samples, metrics match Hail output exactly.

---

### Commit 2 — Add a new WDL task to run coverage DB builder

* Add new WDL task `build_coverage_db` in a new WDL file (do not replace existing yet).

**Acceptance**: WDL runs and produces `coverage.h5` + `coverage_summary.tsv`.

---

### Commit 3 — Shard helper + shard combine task

* Add task to split processed TSV into shard TSVs (by sample rows).
* Add `build_vcf_shard_mt.py` that ingests one shard TSV (containing `gs://` VCF paths) and writes a shard MT tarball.

**Acceptance**: shard MT tarballs produced; covdb hom-ref/DP logic validated vs old pipeline for a 1k cohort.

---

### Commit 4 — Merge tree + finalize-once

* Add merge-group planning (fan-in) and multi-round merges of shard MT tarballs.
* Add `finalize_mt_with_covdb.py` that runs covdb hom-ref/DP imputation and artifact-prone-sites filtering exactly once.

**Acceptance**: full v2 Step 0–3 runs on 50k without OOM; emits a single finalized MT tarball.

---

### Commit 5 — Scale tuning + robustness

* Tune shard size, batch sizes, compression/chunking.
* Add resumability (skip completed shards).

**Acceptance**: stable run on 535k (or staged rollout) with predictable runtime.
