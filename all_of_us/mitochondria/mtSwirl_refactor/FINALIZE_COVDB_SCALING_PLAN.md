# Finalize covdb scaling plan (535k samples)

## Purpose
This document summarizes the current behavior of `determine_hom_refs_from_covdb`, the scaling bottlenecks observed at very large sample counts, and the recommended changes that make the step tractable for ~535k samples while preserving semantics.

## Context
The finalize step updates a sparse mtDNA MatrixTable by using per‑sample coverage from `coverage.h5` to distinguish between **missing** and **hom‑ref** genotypes. The inputs are extremely sparse (observed missing fraction ~0.9993), which means any strategy that scans the entire entry space repeatedly or attempts per‑entry random lookups will not scale.

---

## What the current code does (summary)
The current `determine_hom_refs_from_covdb` implementation:

1. Maps each MT sample to a row index in `coverage.h5` and each MT position to a column index in `coverage.h5`.
2. Splits positions into blocks (`position_block_size`).
3. For each block:
   - Reads coverage for **all samples** × **positions in the block** from HDF5.
   - Builds a per‑block literal and annotates entries with `__cov` using:
     - `mt = mt.annotate_entries(__cov=hl.if_else(mt.__block == block_id, cov_expr, mt.__cov))`
4. After all blocks, applies the hom‑ref logic:
   - If `HL` missing and `DP` > threshold → set `HL=0`, `FT=PASS`, `DP=coverage`.
   - Otherwise keep missing.

### Why this becomes slow
The per‑block `annotate_entries` updates are evaluated **across all entries** each time, even though only one block should change. With ~16–20 blocks, this effectively multiplies entry‑level work ~16–20×. This behavior is the main scaling problem at 535k samples.

---

## Recommended plan (scalable approach)
### A) Keep the same semantics, but avoid repeated full‑MT entry scans
**Key change:** apply entry annotations only to the rows in the current block, then recombine.

**High‑level pattern:**

- Compute `__block` once per row.
- For each block:
  - `mt_b = mt.filter_rows(mt.__block == block_id)`
  - Compute coverage for that block and annotate **only** `mt_b` entries.
  - Apply hom‑ref logic on `mt_b`.
  - Checkpoint `mt_b`.
- Union blocks via a **small fan‑in tree** (not a long linear chain).

This preserves sparsity and ensures each row/entry is processed **once**, not once per block.

### B) Remove unnecessary shuffles
The current `mt.repartition(n_blocks, shuffle=True)` forces a global shuffle. In the block‑local pattern it is not needed and should be removed or replaced with a non‑shuffle repartition only when required for output sizing.

### C) Avoid global `__cov` entry field
Only create `__cov` inside each block MT (`mt_b`). This prevents inflating the entry schema for the full MT and reduces IR size.

### D) Combine with **sample sharding**
At 535k samples, even a perfect block‑local refactor can still be heavy. The strongest scaling improvement comes from **sharding by samples** and running finalize per shard, then `union_cols` the shards and compute cohort statistics afterwards.

Recommended structure:

1. **Split columns** into N shards (e.g., 20–24).
2. For each shard:
   - Run block‑local hom‑ref imputation (A–C above).
3. **Union columns** across shards.
4. Run cohort‑wide row annotations (AC/AN/AF, histograms, hap/pop splits) once.

---

## Why this scales to 535k samples
1. **No repeated full‑matrix entry scans**: each row/entry is processed once.
2. **Sparse‑preserving**: we never densify the MT; missing entries stay missing unless coverage supports hom‑ref.
3. **Reduced driver pressure**: block literals stay small and per‑shard.
4. **Parallelism**: sample shards scale horizontally; each shard’s workload is smaller and independent.
5. **Correctness preserved**: the hom‑ref logic is unchanged; only the execution plan is optimized.

---

## Before vs. After (behavioral differences)
| Aspect | Current code | Recommended code |
|---|---|---|
| Entry‑level evaluation | Full MT per block | Only rows in block |
| Shuffling | Global shuffle (`shuffle=True`) | Removed or minimized |
| `__cov` creation | Whole MT | Block‑local only |
| Recombination | Single MT updated per block | Union of block MTs (fan‑in) |
| Scaling | Work ~ (#blocks × all entries) | Work ~ (all entries once) |

---

## Final notes
- The hom‑ref logic is **identical** to current behavior.
- The plan avoids any “lazy per‑entry lookup” that would cause billions of random HDF5 reads.
- The strategy remains compatible with downstream `add_annotations` and other cohort‑wide statistics, which should be computed **after** unioning shards.

If needed, a follow‑on doc can include concrete code changes or a WDL wiring plan for shard/fan‑in orchestration.
