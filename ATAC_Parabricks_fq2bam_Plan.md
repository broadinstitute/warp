# Plan: Switch ATAC alignment to NVIDIA Parabricks fq2bam (Terra / GPU)

**Repository:** `broadinstitute/warp`
**Branch:** `jw_atac-parabricks-fq2bam` (created off `develop`)
**Status:** Planning complete — no pipeline code changed yet. Awaiting go-ahead to implement.
**Date:** 2026-09-11

> This is a planning artifact for the branch. Delete it before opening the PR (or keep it out
> of the release commit) — it ships no pipeline artifact.

---

## Decisions — LOCKED
1. **GPU type:** T4 (`nvidia-tesla-t4`), driver `535.104.05`, Parabricks `4.5.0-1`,
   `--low-memory` always on (mirrors the scANVI pipeline's proven Terra config).
2. **Strategy:** Augment — add an fq2bam GPU task, **keep** the existing bwa-mem2 CPU task, and
   route by GPU count. **CPU stays the default on every path.** Expose the Parabricks option not
   only in ATAC but through **Multiome and PairedTag** — each caller threads the GPU input into
   its `atac.ATAC` subworkflow call. Additive (existing inputs unchanged), so no contract break.
   (AGENTS.md two-task GPU/CPU pattern; scANVI precedent.)

> **Scope note (Azure):** preserving Azure is **no longer a design driver** — operator confirms
> Azure can be broken and is not needed. *Removing* Azure is out of scope for this work; we just
> stop treating Azure-compat as a constraint. GPU path is GCP-only (`cloud_provider=="gcp"`
> guard); Azure keeps running the CPU path untouched for now.
3. **GPU scaling:** aggregate T4 budget scaled by job size, **hard cap 24 T4**; ≤4 T4/VM (GCP
   limit) realized via a scatter of fq2bam tasks + `MergeSortBam`; host cpu/mem co-scaled;
   `atac_gpu_count` overrides the budget; `gpuType` swappable (A100 path left open).

---

## 0. Straight answer: "will this change our outputs?"
**No, I cannot pre-confirm bit-identical outputs — and it will not be bit-for-bit identical.**
This is an aligner substitution: today's ATAC uses **bwa-mem2 2.2.1** (Open-Omics distributed
build); fq2bam is NVIDIA's accelerated **bwa-mem 0.7.x**. Read placement is >99% concordant,
but ties/edge cases, unmapped-read ordering, the `PA` aux tag, and duplicate marking differ.
The only definitive answer is an empirical same-data diff (see §6, the concordance gate).

Encouragingly, ATAC's own verification **already tolerates aligner nondeterminism**, which is
evidence the team expects run-to-run alignment variation:

| Output | Comparator | Tolerance |
| --- | --- | --- |
| Aligned BAM | `CompareBams` | `mappings_diff_threshold = 0.05` (5%) + lenient header; skip if size differs >200 MB |
| Metrics h5ad | `CompareH5adFilesATAC` | passes at >0.990 similarity |
| Fragment file | `CompareTabix` | md5, else line-count within `max(100, 1e-4 × lines)` |
| Library metrics | `CompareLibraryFiles` | sorted md5, else line diff excluding a few fields (strict-ish) |

So the BAM and h5ad have real headroom; **fragment-file line count and library-metric values**
are the exposed surfaces. If they drift beyond tolerance we either widen those two tolerances
in `verification/**` (allowed — no version bump) or regenerate truth.

---

## 1. How ATAC alignment works today (invariants to preserve)
1. `FastqProcessATAC` (warp-tools `fastqprocess`) corrects the 10x barcode and writes it into
   the **FASTQ comment**: `@read CR:Z:<raw> CB:Z:<corrected> CY:Z:<qual>`
   (warp-tools `tools/fastqpreprocessing/src/fastq_common.cpp:156,185`).
2. `TrimAdapters` (cutadapt) trims per shard (scatter).
3. `BWAPairedEndAlignment` (`pipelines/wdl/atac/atac.wdl:358`) runs distributed **bwa-mem2** with
   `PARAMS="+R '@RG…' +C"`. The `+C` = bwa `-C` copies the FASTQ comment into the BAM, so the
   **`CB` tag** is present. Output is coordinate-sorted, reheadered with an `@CO` reference line.
4. `CreateFragmentFile` (SnapATAC2) reads the BAM via `barcode_tag="CB"` (or `BB` when
   `preindex`) → fragment file + `metrics.h5ad` + `library_metrics.csv`.

**⇒ Non-negotiable invariant:** the aligned BAM must still carry the `CB` tag (and `BB` for the
paired-tag preindex path). fq2bam must forward `-C` to bwa via `--bwa-options "-C"`. Whether
fq2bam's GPU sort/output path preserves comment tags is the **#1 thing to verify on real data**
before anything else — the entire downstream fragment/metrics step depends on it.

---

## 2. Blast radius (why the design is additive)
- ATAC is imported as a sub-workflow by **Multiome** (`pipelines/wdl/multiome/Multiome.wdl:117`)
  and **PairedTag** (`pipelines/wdl/paired_tag/PairedTag.wdl:127`); both thread
  `num_threads_bwa`, `mem_size_bwa`, `cpu_platform_bwa`, `vm_size` into ATAC.
- Those keys are set across ~12 test JSONs in `atac/`, `multiome/`, `paired_tag/`.
- WARP's sub-workflow input contract means removing an ATAC input forces edits to every caller +
  every test JSON. **So keep existing inputs; add GPU inputs with defaults.**
- Editing `atac.wdl` cascades version bumps to **ATAC, Multiome, PairedTag** (each: a
  `pipeline_version` bump + changelog entry). Confirm the exact set with
  `scripts/validate_release.sh -g origin/staging`.
- **Expose the GPU option in callers (required):** add the GPU routing input to `Multiome.wdl`
  and `PairedTag.wdl` and forward it into their `atac.ATAC` call; default unset → CPU. Purely
  additive → no contract break.

---

## 3. Design — additive two-task routing + scalable GPU scatter
Cromwell rejects `gpuCount:0` and optional-typed GPU attributes, so a task is all-GPU or
all-CPU (per AGENTS.md; scANVI `MultiomeLabelTransfer`/`…Cpu` is the precedent).

1. **New GPU task** `BWAPairedEndAlignmentParabricks` in `atac.wdl`:
   - `runtime`: `docker: "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"`, `gpuType`,
     `gpuCount`, `nvidiaDriverVersion: "535.104.05"`, `zones`, plus cpu/memory/disk/bootDiskSizeGb
     sized to Parabricks host minimums.
   - `command`: extract reference; run
     ```
     pbrun fq2bam --ref genome.fa \
       --in-fq R1 R3 "@RG\tID:RG1\tSM:RGSN1" \
       --bwa-options "-C" \
       --low-memory \
       --no-markdups \
       --out-bam <id>.bam
     ```
     (`--bwa-options "-C"` preserves `CB`; `--no-markdups` matches today's no-dup-marking;
     `--low-memory` required on ≤48 GB GPUs). Re-apply the `@CO` reheader for header parity.
   - Keep outputs identical (`bam_aligned_output` + a logs tarball).
2. **Keep** the existing `BWAPairedEndAlignment` (bwa-mem2 CPU) task unchanged.
3. **Route** in the workflow (default = CPU):
   ```
   if (gpu_count > 0)  → scatter BWAPairedEndAlignmentParabricks over shards → MergeSortBam
   if (gpu_count == 0) → BWAPairedEndAlignment (existing bwa-mem2 CPU path)
   select_first([...]) gathers the final BAM
   ```
4. **Validate cloud** — require `cloud_provider == "gcp"` when `gpu_count > 0` via
   `utils.ErrorWithMessage` (GPU path is GCP/Terra-only; Azure keeps the CPU path — Azure is no
   longer a design constraint, but GPU stays GCP-only).
5. **Callers (Multiome, PairedTag)** — add `atac_gpu_count` (+ any GPU knobs) to each caller's
   inputs and forward into the `atac.ATAC` (`Atac` / `Atac_preindex`) call. Default unset = CPU.
   Version bumps already cascade to both. Keep the CPU path selectable everywhere.

---

## 4. GPU scaling — T4, count scaled by job size, aggregate ceiling 24
Parabricks needs ≥16 GB VRAM; 16–48 GB GPUs require `--low-memory` (T4 = 16 GB).

**Hard GCP constraint:** T4 is capped at **4 GPUs per VM** (valid 1/2/4), and one Cromwell task
= one VM. The largest single-VM GPU box on GCP is 16 (A100 `a2-megagpu-16g`). So >4 T4 cannot
live on one node — reaching up to 24 requires **scatter across VMs + merge**.

**Scalable design (scatter fq2bam + merge), user ceiling = 24 T4 aggregate:**
```
# aggregate T4 budget from total trimmed-FASTQ size, hard-capped at 24 (overridable)
Int gpu_budget     = min(24, max(1, ceil(total_fastq_gib / GIB_PER_GPU)))   # GIB_PER_GPU tuned in §6
Int gpu_per_shard  = min(4, gpu_budget)                 # GCP T4 per-VM cap = 4
Int n_align_shards = ceil(gpu_budget / gpu_per_shard)   # e.g. 24 -> 6 shards x 4 T4
# each fq2bam shard task: T4 x gpu_per_shard, --low-memory, --bwa-options "-C", --no-markdups,
#   host co-scaled  cpu ≈ 12 x gpu_per_shard,  mem ≈ 50 x gpu_per_shard GB
#   (Parabricks host mins: 2 GPU -> ≥100 GB / 24 thr, 4 GPU -> ≥196 GB / 32 thr)
# then MergeSortBam over the per-shard coordinate-sorted BAMs -> final BAM (order-only; keeps CB)
```
- `n_align_shards` drives `FastqProcessATAC.num_output_files` on the GPU path.
- Plumbing → 1 T4 / 1 shard / no merge. Large libraries → up to 6 shards × 4 T4, merged.
- `atac_gpu_count` overrides the budget; `gpuType` stays swappable for a future non-T4 run.
- `MergeSortBam` is already imported into `atac.wdl` (currently a dead import per AGENTS.md
  known-debt list) — this puts it to use.
- **Caveat:** 24 concurrent T4 needs matching project/zone **GPU quota**; document it.
- **Alternative (not chosen):** single **A100** VM (≤16/VM, no `--low-memory`, no scatter) —
  fewer, faster, pricier GPUs. Left reachable via overridable `gpuType`/driver.

---

## 5. Reference index
fq2bam uses classic **bwa-mem 0.7.x**, which needs a `.amb/.ann/.bwt/.pac/.sa` index — **not**
the bwa-mem2 index (`.bwt.2bit.64`, `.0123`, …) inside the current
`bwa-mem2-2.2.1-*.tar` reference tars used by the tests. Options:
- (a) Point GPU runs at a classic-bwa tar (Broad publishes `…/hg38/v0/bwa/…`), or
- (b) Let fq2bam build the index from the `.fa` at runtime (adds time).

GPU test cases must use a compatible reference.

---

## 6. Test runtime updates & validation gates
**Test runtimes** ("update all test runtimes"): for GPU-enabled cases, the CPU/BWA knobs
(`num_threads_bwa`, `mem_size_bwa`, `cpu_platform_bwa`, `vm_size`) no longer drive the aligner —
add the GPU inputs (`atac_gpu_count`, etc.) and a classic-bwa `tar_bwa_reference` instead.
- `pipelines/wdl/atac/test_inputs/Plumbing/{10k_pbmc_downsampled,8k_cortex}.json` → `atac_gpu_count = 1`.
- **Multiome + PairedTag:** expose the option; add a GPU case for each **and keep a CPU case**;
  forward the new input through `TestMultiome.wdl` / `TestPairedTag.wdl` and their test JSONs.
- Forward any new top-level input through `verification/test-wdls/TestATAC.wdl`.

**Gates:**
1. `womtool validate` on `atac.wdl` + every importer (Multiome, PairedTag, TestATAC,
   TestMultiome, TestPairedTag) and the TestATAC chain.
2. **Concordance run** on Terra: fq2bam vs bwa-mem2 on `10k_pbmc_downsampled` → diff
   **CB-tag presence**, BAM mapping % / flagstat, fragment line counts, h5ad + library metrics.
   This answers §0 with numbers and tunes `GIB_PER_GPU`.
3. If within tolerance → keep truth (optionally widen fragment/library tolerances in
   `verification/**` — no version bump). Else regenerate truth via
   `workflow_dispatch updateTruth:true` on the right `truthBranch`.
4. `scripts/validate_release.sh -g origin/staging -i true` until all importers are green.
5. Changelogs: ATAC (+ Multiome, PairedTag) bumped with dated entries.
6. Docs: update the ATAC docs-site page (new GPU inputs); `yarn --cwd=website build`.

---

## 7. Deliverables checklist (post-approval)
- [ ] `atac.wdl`: new Parabricks GPU task + scatter/merge + `gpu_count`-based routing + scaled
      GPU/host inputs; keep the bwa-mem2 CPU path; `cloud_provider=="gcp"` guard for GPU.
- [ ] `Multiome.wdl` + `PairedTag.wdl`: add GPU input, forward into the ATAC call; CPU default.
- [ ] `pipeline_version` bump + changelog entry: ATAC, Multiome, PairedTag.
- [ ] Test JSON runtime updates (+ classic-bwa reference for GPU cases), incl. a GPU + a CPU
      case for Multiome and PairedTag.
- [ ] `TestATAC.wdl` / `TestMultiome.wdl` / `TestPairedTag.wdl` forward any new top-level inputs;
      CI `paths:` still correct.
- [ ] womtool + `validate_release.sh` pass; docs build.
- [ ] Concordance report; truth updated only if justified.

---

## Key file references
- Pipeline: `pipelines/wdl/atac/atac.wdl` (workflow + `BWAPairedEndAlignment` at line 358)
- Barcode source: warp-tools `tools/fastqpreprocessing/src/fastq_common.cpp:156,185`
- GPU pattern precedent: `pipelines/wdl/scanvi/scANVI.wdl:470` (`MultiomeLabelTransfer`/`…Cpu`)
- Callers: `pipelines/wdl/multiome/Multiome.wdl:117`, `pipelines/wdl/paired_tag/PairedTag.wdl:127`
- Verification: `verification/VerifyATAC.wdl`, `verification/VerifyTasks.wdl`
- Test wrapper / CI: `verification/test-wdls/TestATAC.wdl`, `.github/workflows/test_atac.yml`
- Conventions: `AGENTS.md` (GPU two-task pattern, cascading version bumps, test registration)

## External references
- fq2bam tool reference — https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam
- Parabricks install requirements (GPU memory, `--low-memory`, drivers) —
  https://docs.nvidia.com/clara/parabricks/get-started/installation-requirements
- Parabricks on Terra — https://docs.nvidia.com/clara/parabricks/tutorials/cloud-usage-guides/terra
- Official Parabricks WDLs (fq2bam) — https://github.com/clara-parabricks-workflows/parabricks-wdl
- GCP GPUs (T4 max 4/VM) — https://docs.cloud.google.com/compute/docs/gpus/about-gpus
- GCP accelerator-optimized machines (A100 up to 16/VM) —
  https://docs.cloud.google.com/compute/docs/accelerator-optimized-machines
