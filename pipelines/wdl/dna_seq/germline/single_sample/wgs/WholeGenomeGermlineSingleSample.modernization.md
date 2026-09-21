# Modernization & Optimization Assessment — WholeGenomeGermlineSingleSample (WARP)

**Target:** `pipelines/wdl/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl` (pipeline_version 3.3.8) and its imported task graph.
**Date:** 2026-09-21
**Author:** Claude (analysis for jway@broadinstitute.org)
**Scope of this document:** findings + feasibility, *not* a hand-off implementation plan. Ground rules applied throughout: **inputs unchanged**, **outputs scientifically equivalent (not byte-identical)** unless published data shows a new approach is measurably better.

---

## 1. What the pipeline is today

The top-level WDL wires together five imported subworkflows/task libraries:

```
WholeGenomeGermlineSingleSample.wdl
├── UnmappedBamToAlignedBam.wdl        # align + dedup + sort + BQSR
│     ├── Alignment.wdl                # SamToFastq | bwa mem | MergeBamAlignment  (streamed)
│     ├── DragmapAlignment.wdl         # alt aligner (DRAGMAP / dragen-os)
│     ├── SplitLargeReadGroup.wdl      # split RGs > 20 GiB, align shards, regather
│     ├── BamProcessing.wdl            # SortSam, MarkDuplicates, BaseRecalibrator, ApplyBQSR, Gather*, CheckContamination
│     └── Qc.wdl                       # yield / readgroup QC
├── AggregatedBamQC.wdl                # Picard multi-metrics, fingerprinting
├── Qc.wdl                             # CollectWgsMetrics / CollectRawWgsMetrics
├── BamToCram.wdl                      # final CRAM
└── VariantCalling.wdl                 # HaplotypeCaller (GATK3 *or* GATK4/DRAGEN), Reblock, MergeVcfs, ValidateVCF, metrics
```

Data flow on the default path: per-read-group uBAM → **BWA-MEM** (streamed through Picard SamToFastq → bwa → MergeBamAlignment) → aggregate **Picard MarkDuplicates** (queryname-assumed) → **Picard SortSam** (coordinate) → **GATK BQSR** (scattered BaseRecalibrator + ApplyBQSR) → **HaplotypeCaller** (scattered) → MergeVcfs → ReblockGVCF → validate/metrics. Final aligned output is CRAM; final variant output is a reblocked GVCF.

### 1.1 Tool/version inventory (core path)

| Step | Tool | Version in use | Container | Age |
|---|---|---|---|---|
| Alignment | BWA-MEM | **0.7.15** | `samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748` | 2016 |
| SamToFastq / MergeBamAlignment / SortSam / MarkDuplicates / Gather / QC metrics | Picard | **2.26.10** | `picard-cloud:2.26.10` | Dec 2021 |
| BQSR (BaseRecalibrator / ApplyBQSR / GatherBQSRReports) | GATK | **4.3.0.0** | `broad-gatk/gatk:4.3.0.0` | 2022 |
| HaplotypeCaller (GATK4 path) / Reblock / ValidateVCF / DRAGSTR | GATK | **4.6.1.0** | `broad-gatk/gatk:4.6.1.0` | 2024 |
| HaplotypeCaller (**default** path) | **GATK 3.5** (`GATK35.jar`) wrapped by GATK4 PrintReads | `gatk:1.3.0-4.2.6.1-1649964384` | 3.5 ≈ 2015 |
| Alt aligner | DRAGMAP / dragen-os | 1.2.1 | `dragmap:1.1.2-1.2.1-2.26.10-...` | opt-in |
| Contamination | VerifyBamID (verifyBamID2) | 1.0.1 | `verify-bam-id:1.0.1-...-1639071840` | 2019 |
| bamout merge | **samtools 1.3.1** | 1.3.1 | `biocontainers/samtools:1.3.1` | 2016 |
| contamination subsetting | bedtools | 2.27.1 | `bedtools:2.27.1` | 2017 |

**These BAM-processing tasks are shared with `ExomeGermlineSingleSample.wdl`** — any change to `BamProcessing.wdl` / `UnmappedBamToAlignedBam.wdl` affects both pipelines. That doubles the payoff and the blast radius.

---

## 2. Headline findings

Two things stand out before the itemized list:

1. **The default variant caller is GATK 3.5.** At the top level, `use_gatk3_haplotype_caller = true` (line 59). That routes to `HaplotypeCaller_GATK35_GVCF`, which literally runs `java -jar /usr/gitc/GATK35.jar -T HaplotypeCaller … --max_alternate_alleles 3 -variant_index_type LINEAR --read_filter OverclippedRead` after a GATK4 `PrintReads` shim. A ~2015 caller is the out-of-the-box behavior in a 2026 pipeline. The modern GATK4 caller already exists in the same file and is used when the flag is flipped. This is the single largest "science is not as good as it could be" item, and it is almost free to change.
2. **The pipeline already contains its own modernization on-ramps** — a GATK4 path, a DRAGEN-GATK path (`--dragen-mode`, DRAGSTR auto-calibration, DRAGEN hard-filtering), and a DRAGMAP aligner path — all gated behind default-off flags. And the organization has **already built the reference assets and the engineering patterns** for the bigger swaps (bwa-mem2 indices in the public bucket; a Parabricks fq2bam GPU task in ATAC/Multiome/PairedTag). Much of what follows is "promote what you already have / reuse what you already built," not greenfield.

---

## 3. Efficiency & performance opportunities (output stays scientifically equivalent)

### 3.1 Replace BWA-MEM 0.7.15 with BWA-MEM2 — *drop-in, index already exists*
**Now:** `bwa mem -K 100000000 -p -v 3 -t 16 -Y` on BWA **0.7.15** (`Alignment.wdl:52`, `bwa_commandline` at `UnmappedBamToAlignedBam.wdl:55`), 16 vCPU / 14 GiB.
**Why suboptimal:** BWA-MEM2 is the maintained successor; it is a functionally-equivalent re-engineering of the same algorithm delivering roughly **1.5–2×** wall-clock on the alignment step for identical output. 0.7.15 is also two point releases behind the last classic BWA (0.7.17/0.7.18).
**Improvement / feasibility:**
- BWA-MEM2 2.2.1 is documented as producing output **identical to BWA-MEM 0.7.17** and is a drop-in replacement (~50–100% seeding speedup) ([Vasimuddin et al., IPDPS 2019](https://ieeexplore.ieee.org/document/8820962)).
- **The hg38 bwa-mem2 index already exists in Broad's own public bucket:** `gs://gcp-public-data--broad-references/hg38/v0/bwa/v2_2_1/bwa-mem2-2.2.1-Human-GENCODE-build-GRCh38.tar` (used today by Multiome/PairedTag/ATAC). No new reference-generation project needed for the human default.
- Keep the same `-K 100000000` chunk size to preserve determinism; the streamed SamToFastq→aligner→MergeBamAlignment structure is unchanged.
- **Caveat to validate:** bwa-mem2 matches bwa-mem **0.7.17**, and this pipeline is pinned to **0.7.15**. Expect near-identical but not historically-identical placement vs the current 0.7.15 truth set — so this needs a concordance gate, not a blind swap. Memory footprint of the bwa-mem2 index is larger (index ~2× on disk); the 14 GiB task will need a bump.

### 3.2 Optional GPU alignment via NVIDIA Parabricks fq2bam — *precedent already in this repo*
**Why relevant:** The team has **already implemented** a Parabricks `fq2bam` GPU alignment task in `atac.wdl` (`clara-parabricks:4.5.0-1`, T4, scatter ≤4 GPU/VM + merge, GCP-only, additive/default-off) and exposed it through Multiome and PairedTag. The design reasoning in `ATAC_Parabricks_fq2bam_Plan.md` transfers almost verbatim to WGS.
**Improvement / feasibility:**
- fq2bam = accelerated BWA-MEM + coordinate sort (+ optional MarkDuplicates + BQSR) that follows GATK Best Practices. NVIDIA's GIAB benchmarking reports **>0.9999 precision/recall equivalence** to the BWA+MarkDuplicates+BQSR+HaplotypeCaller steps of *this very WARP workflow*; an independent AWS/HG001 30× study reports **hap.py F1 = 0.9999** vs GATK4 ([NVIDIA fq2bam docs](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam.html), [AWS HPC blog](https://aws.amazon.com/blogs/hpc/benchmarking-the-nvidia-clara-parabricks-germline-pipeline-on-aws/)).
- A 30× WGS BAM in minutes vs hours; the bigger win is that fq2bam can *fuse* align+sort+markdup+BQSR, collapsing four of this pipeline's steps into one task (see 3.3/3.4).
- Well-documented, controllable sources of non-identity: pass matching `--bwa-options="-K 100000000"`; the `PA` tag ordering/rounding and unmapped-read sort order differ. Use the same two-task GPU/CPU additive pattern (CPU stays default) the team already adopted — no input-contract change.
- **Feasibility notes:** GCP-only (`gpuType`/`gpuCount` runtime attrs), license/registry pull for `nvcr.io`, and cost model differs (fewer, pricier VM-hours). Best framed as an opt-in fast lane, exactly like ATAC.

### 3.3 Picard SortSam is single-threaded — the clearest pure-speed win
**Now:** `BamProcessing.SortSam` runs `picard SortSam … MAX_RECORDS_IN_RAM=300000` on **1 CPU** (`BamProcessing.wdl:19-62`). The task's own comment notes it "spills to disk a lot." It runs once over the *entire aggregated sample* BAM.
**Why suboptimal:** Picard's sort core is single-threaded; adding RAM doesn't help (documented behavior). `samtools sort -@ N` scales across cores (≈3–4× realistic on many cores) and is the standard fast coordinate sorter ([biowdl ngs-performance-choices](https://github.com/biowdl/ngs-performance-choices), [Picard #529](https://github.com/broadinstitute/picard/issues/529)).
**Improvement / feasibility:**
- This SortSam is *downstream* of MarkDuplicates and only needs **coordinate** order — precisely the case where `samtools sort` is a clean substitute (the queryname-sort tie-break differences between the tools do **not** apply to coordinate sort). Output is scientifically equivalent (same coordinate order; index + md5 regenerated with `samtools index` / `md5sum`).
- Feasible as an internal-tool swap with no interface change. Multi-thread the task (e.g. 8 vCPU) and it stops being a serial tentpole on the critical path.

### 3.4 Collapse sort + mark-duplicates (+ BQSR) into one pass — elprep
**Now:** MarkDuplicates (single-thread, `BamProcessing.wdl:66`), then SortSam (single-thread), then scattered BQSR are three separate tasks each re-reading/re-writing the whole-sample BAM.
**Why suboptimal:** Repeated whole-genome I/O and JVM spin-up per step; MarkDuplicates is also single-threaded.
**Improvement / feasibility:**
- **elprep** runs mark-duplicates + sort + BQSR in a single (in-memory or `sfm`) pass and is documented to produce **the same BAM, metrics, and recalibration files as GATK4**, at ~**7×** the throughput in sfm mode ([elprep docs/paper](https://github.com/ExaScience/elprep)). This is the "functionally equivalent by design" option for the whole pre-processing block.
- Feasibility: elprep equivalence claims are for the GATK4 flavors of these steps; it needs a memory-heavy (or sfm disk-backed) VM and its own concordance validation, but it is the highest-leverage single change for pre-processing cost. Lower-risk partial version: just multi-thread MarkDuplicates' successor (samtools markdup / sambamba) — but duplicate-flag *sets* can differ between markers, so that path needs the concordance gate; elprep is the option with the strongest published equivalence.

### 3.5 Container & version drift / inconsistency — cheap hygiene
**Findings:**
- **Picard 2.26.10 (Dec 2021)** is used in ~15 places; current Picard is 3.x. Several CVE/JDK and bug fixes since.
- **Two GATK versions on one path:** BQSR pinned to **4.3.0.0** while HaplotypeCaller/Reblock use **4.6.1.0**. No reason for BQSR to lag; unify on 4.6.x (BQSR tables are stable across these minor versions — low equivalence risk, easy to validate).
- **samtools 1.3.1 (2016)** in `VariantCalling.MergeBamouts` and **1.11** in `AggregatedBamQC` — pin both to a current samtools. (MergeBamouts only runs when `make_bamout=true`, so low urgency but it's a 10-year-old container.)
- bedtools 2.27.1, VerifyBamID 1.0.1 — review for currency.
**Why suboptimal:** security surface, reproducibility ambiguity, and losing free performance/bugfix improvements. **Improvement:** a coordinated container-bump pass with per-tool concordance checks. Individually low-risk; collectively meaningful.

#### 3.5.1 Picard 2.26.10 → 3.x specifically — mechanically easy; the Java jump is already solved
The Picard bump is the most-asked and the least scary. Picard 3.x requires **Java 17** (3.0.0 dropped Java 8), which sounds like a blocker but is not: the JRE is bundled in the container and the WDLs invoke a generic `java -jar`, so the Java version is the image's concern, not the WDL's. **The org already builds and runs `us.gcr.io/broad-gotc-prod/picard-cloud:3.0.0`** (in `RNAWithUMIsTasks.wdl`, `atac.wdl`, `verification/VerifyMetrics.wdl`) — a proven Java-17 Picard 3.0.0 image, at the **same jar path** (`/usr/picard/picard.jar`). Picard 3.x also still defaults to the legacy `KEY=VALUE` argument parser, so existing command lines (`SortSam INPUT=… SORT_ORDER=coordinate`, `-Dsamjdk.compression_level`, `-Xms/-Xmx`) run **unchanged**.

Difficulty splits three ways:

| Tier | Tasks (WGS core path) | Change required | Effort |
|---|---|---|---|
| **Trivial** | ~46 calls pinned to `picard-cloud:2.26.10` — SortSam, MarkDuplicates, Gather\*, all `Qc.wdl`/`Metrics.wdl` metrics, MergeVcfs | swap tag string `2.26.10` → `3.0.0` | string edit; image already exists |
| **Docker rebuild** | ~11 calls on `/usr/gitc/picard.jar` — **`Alignment.wdl`** (SamToFastq + MergeBamAlignment + SamSplitter), `Utilities.wdl`, `UltimaGenomicsWholeGenomeGermlineTasks.wdl` | Picard is fused into the `samtools-picard-bwa` combo image; bumping it means rebuilding that combo container | real work; on the critical path |
| **Minor / opt-in** | 1 call on `/picard/picard.jar` (Dragmap combo image) | same combo-rebuild story | low urgency |

Only the combo-image tier is an actual engineering lift — you cannot retag it; the Picard jar ships inside a bwa+picard+samtools image. A safe first move is to bump only the `picard-cloud` tier (all the sort/dedup/metrics tasks) and leave `Alignment.wdl` on 2.26.10 — its Picard steps only merge the uBAM; coordinate sorting, dedup and metrics all happen downstream on `picard-cloud`.

**The real cost is output equivalence, not the edit.** 2.26.10 → 3.0.0 is a large **htsjdk** jump. Concordance-gate the tools whose *values* can move:
- **MarkDuplicates** — duplicate-flag set and library-size metrics can shift; highest risk against the "scientifically equivalent" bar.
- **CollectWgsMetrics / CollectQualityYieldMetrics / insert-size / GC-bias** — numeric QC metrics may differ at the margins (lower stakes; these are QC outputs, not the callset).
- **SortSam / GatherBamFiles** — order/concatenation only → safe.

Net: an afternoon of tag edits for the picard-cloud tier plus a GIAB before/after diff on MarkDuplicates metrics and the QC metric files; the combo-image rebuild is a separate, larger chunk taken only if newer SamToFastq/MergeBamAlignment behavior is actually wanted. (A `picard-cloud:2.26.11` image also already exists as a smaller interim step if a full 3.x jump is deemed too much validation at once.)

### 3.6 Runtime knobs — smaller, lower-confidence tunings
- **HDD everywhere.** Most tasks request `local-disk … HDD`. The I/O-bound serial steps (SortSam, MarkDuplicates, Gather) are exactly where SSD would cut wall-clock; worth A/B-ing cost vs time.
- **`compression_level = 2`** for intermediates (`UnmappedBamToAlignedBam.wdl:57`) is a deliberate speed-over-size choice and is probably fine; revisit only alongside a sort/dedup rework.
- **Scatter granularity:** BQSR and HaplotypeCaller scatter counts are static (`bqsr_divisor`, `hc_divisor` heuristics). Not wrong, but a candidate for right-sizing once the aligner/pre-processing costs move.
- These are "confirm with a cost experiment" items, not clear wins — listed for completeness.

---

## 4. Scientific opportunities (output changes; each backed by published data)

### 4.1 Stop defaulting to GATK 3.5 HaplotypeCaller — *the biggest science gap, near-free to close*
**Now:** default `use_gatk3_haplotype_caller = true` → `GATK35.jar -T HaplotypeCaller`.
**Why suboptimal:** GATK 3.5 predates a decade of assembly/genotyping model improvements. The GATK4 caller is present in the same repo and is what every current Best-Practices recommendation uses. The only reason to keep 3.5 is historical byte-continuity with legacy callsets.
**Improvement / feasibility:**
- Flip the default to the GATK4 path (`use_gatk3_haplotype_caller = false`). The GATK4 vs GATK3 concordance is high on easy sites but GATK4 is materially better on indels and hard regions; the "functional equivalence" framework ([Regier et al., *Nat Commun* 2018](https://www.nature.com/articles/s41467-018-06159-4)) exists precisely so that upstream harmonization lets you innovate on the caller without breaking joint analysis.
- **Feasibility:** trivial mechanically (flag already threaded end-to-end; a guard already forbids GATK3 + DRAGEN together). The real work is the callset-migration decision + a GIAB concordance report, not code. This should be recommendation #1.

### 4.2 Promote DRAGEN-GATK (functional-equivalence / max-quality modes) from opt-in to first-class
**Now:** `dragen_functional_equivalence_mode` / `dragen_maximum_quality_mode` exist and orchestrate DRAGMAP alignment + `--dragen-mode` HaplotypeCaller + DRAGSTR + DRAGEN hard-filtering, but default off; classic BWA+BQSR+GATK is the default.
**Why suboptimal:** DRAGEN-GATK's STR-aware genotyping (DRAGSTR) improves indel/STR accuracy and lets you **drop BQSR entirely** (the modes already set `perform_bqsr=false`), removing two scatter stages. It is the direction GATK itself has taken.
**Improvement / feasibility:** the machinery is built and validated enough to ship behind a flag; the question is whether to make it (or the max-quality preset) the recommended default. Needs a head-to-head GIAB benchmark on Broad data and a callset-continuity decision. Feasible now — it's a defaulting/validation decision, not new engineering.

### 4.3 Consider DeepVariant as an alternative caller for accuracy-first use cases
**Why:** Across independent benchmarks and PrecisionFDA Truth Challenges, DeepVariant matches GATK on SNPs at high coverage (F1 ≈ 0.98–0.99) and **beats HaplotypeCaller on indels, especially in low-complexity/homopolymer regions**; at cohort scale DeepVariant+GLnexus showed markedly lower Mendelian-violation and F1-error rates ([Yun et al., *Bioinformatics* 2020](https://academic.oup.com/bioinformatics/article/36/24/5582/6064144); [Supernat et al., *Sci Rep* 2018](https://www.nature.com/articles/s41598-018-36177-7); [Lin et al., *Sci Rep* 2022](https://www.nature.com/articles/s41598-022-05833-4)).
**Feasibility / caveats:** This is the largest scientific departure. DeepVariant produces its own GVCF (cohort merging is via GLnexus, not GenomicsDB/GATK joint genotyping), so it is **not** functionally equivalent to the current GVCF and would branch the downstream joint-genotyping story. Best positioned as an *alternative pipeline / opt-in caller* for accuracy-first projects, explicitly justified by the published indel advantage — not a swap under the "scientifically equivalent" rule. Parabricks also ships a GPU DeepVariant, so it composes with 3.2.

### 4.4 Smaller science-adjacent items
- **VerifyBamID (contamination):** 1.0.1 (2019). Confirm it's the latest verifyBamID2 and that the contamination-resource panels are current; contamination feeds directly into HaplotypeCaller's `-contamination`, so drift here is scientific, not cosmetic.
- **`--read_filter OverclippedRead` / `--max_alternate_alleles 3`** on the GATK3 path are legacy defaults; they disappear naturally when 4.1 lands.

---

## 5. Feasibility summary

| # | Change | Output impact | Published/precedent evidence | Main blocker to confirm |
|---|---|---|---|---|
| 3.1 | BWA-MEM → **BWA-MEM2** | Near-identical (matches bwa 0.7.17, not 0.7.15) | IPDPS 2019; identical-output claim | hg38 index exists ✓; validate vs 0.7.15 truth; +RAM/disk |
| 3.2 | Optional **Parabricks fq2bam** GPU | >0.9999 concordant | NVIDIA GIAB + AWS HG001; **in-repo ATAC precedent** ✓ | GPU/GCP-only, cost, `-K` match |
| 3.3 | **samtools sort** for coordinate SortSam | Equivalent (coordinate order) | biowdl bench; Picard #529 | trivial; regen index/md5 |
| 3.4 | **elprep** fused sort+markdup(+BQSR) | Same BAM/metrics/recal as GATK4 | elprep 7.4× w/ identical files | RAM-heavy VM; concordance gate |
| 3.5 | Container/version unification | Equivalent | — | coordinated bump + checks |
| 3.6 | SSD / scatter / compression tuning | Equivalent | — | cost experiment |
| 4.1 | **Default → GATK4** HaplotypeCaller | Better indels/hard regions | Regier 2018 FE framework | callset migration decision |
| 4.2 | Promote **DRAGEN-GATK** default | Better STR/indel; drops BQSR | built-in; GATK direction | GIAB head-to-head + continuity |
| 4.3 | **DeepVariant** (opt-in) | Better indels; **not** FE-equivalent | pFDA, Sci Rep, Bioinformatics | branches joint-genotyping |
| 4.4 | VerifyBamID / resources currency | Potentially scientific | — | version/panel audit |

---

## 6. Explicitly out of scope / not recommended (and why)

- **Changing inputs** (uBAM contract, reference layout) — excluded by the ground rules. All recommendations above keep the input contract.
- **Ripping out BQSR unconditionally** — only justified on the DRAGEN-GATK path (4.2), which already handles it. For the classic BWA+GATK4 path, BQSR stays.
- **Byte-identical guarantees** — none of the aligner/caller swaps are byte-identical; every one needs a concordance gate. The pipeline's sibling verification harness (and ATAC's existing tolerance-based comparators) is the right place to set those gates.
- **A step-by-step implementation plan** — deliberately omitted per the request. The next artifact would be: pick the 2–3 changes to pursue, then write per-change concordance-test plans (GIAB HG001/HG002 truth, hap.py F1 gates) and container-build tasks.

---

## 7. Suggested sequencing (if pursued)

Cheapest-highest-value first: **4.1 (GATK4 default)** and **3.3 (samtools sort)** and **3.5 (version unification)** are low-risk, high-clarity. **3.1 (bwa-mem2)** is next given the index already exists. **3.4 (elprep)** and **4.2 (DRAGEN-GATK default)** are the structural bets that need real benchmarking. **3.2 (Parabricks)** and **4.3 (DeepVariant)** are opt-in lanes that reuse work/tools already proven elsewhere in the org.

---

## Sources
- Vasimuddin et al., "Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems," IPDPS 2019 — https://ieeexplore.ieee.org/document/8820962 ; BWA-MEM2 repo — https://github.com/bwa-mem2/bwa-mem2
- NVIDIA Parabricks fq2bam docs — https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam.html ; AWS HPC benchmarking blog — https://aws.amazon.com/blogs/hpc/benchmarking-the-nvidia-clara-parabricks-germline-pipeline-on-aws/
- Regier et al., "Functional equivalence of genome sequencing analysis pipelines…," Nature Communications 2018 — https://www.nature.com/articles/s41467-018-06159-4
- Yun et al., "Accurate, scalable cohort variant calls using DeepVariant and GLnexus," Bioinformatics 2020 — https://academic.oup.com/bioinformatics/article/36/24/5582/6064144
- Supernat et al., "Comparison of three variant callers for human whole genome sequencing," Sci Rep 2018 — https://www.nature.com/articles/s41598-018-36177-7
- Lin et al., "Comparison of GATK and DeepVariant by trio sequencing," Sci Rep 2022 — https://www.nature.com/articles/s41598-022-05833-4
- elprep — https://github.com/ExaScience/elprep ; biowdl ngs-performance-choices — https://github.com/biowdl/ngs-performance-choices ; Picard SortSam speed issue #529 — https://github.com/broadinstitute/picard/issues/529
- In-repo precedent: `ATAC_Parabricks_fq2bam_Plan.md`; `pipelines/wdl/atac/atac.wdl` (BWAPairedEndAlignmentParabricks); bwa-mem2 index `gs://gcp-public-data--broad-references/hg38/v0/bwa/v2_2_1/bwa-mem2-2.2.1-Human-GENCODE-build-GRCh38.tar`
