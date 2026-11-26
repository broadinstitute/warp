# AoU + AnVIL chrX panel: running everything end to end

Status of the code in this bundle (24 November 2026):

| Piece | Status |
|---|---|
| Steps 1 and 2 (`Step1_2_Site_Stats.py`) | run on real data, 200 kb and 5 Mb; output checked |
| `Eval_Callset.py` and normalisation | run on real data, 5 Mb; output not yet used |
| `check_nonpar_representation.py`, `infer_sex_from_vds.py` | run on real data, both cohorts |
| Step 3 (`Step3_Filter_Rename_Combine.py`) | first version ran but stalled in the combiner; the fixed version is untested |
| Step 4 (`Step4_Densify_Export.py`), masking, AoU-only arm | never run |
| `filter_sweep.py` | tested on synthetic data only |
| WDL | runs in Terra for Steps 1, 2 and the evaluation callset; not validated with `womtool` |

Finish the smoke test in section 2 (Steps 3 and 4) before the 5 Mb panel run.

Timings so far on 2 + 8 `n2d-highmem-8` workers, 5 Mb: Step 1 45 min, Step 2 20 min (in parallel),
evaluation callset 20 min.

## 0. What is in the bundle

| File | Runs where | Purpose |
|---|---|---|
| `pipeline/hail_in_wdl_optimized.wdl` | Terra | Hail workflow (replaces `hail_in_wdl_first_try.21.wdl`) |
| `pipeline/Step1_2_Site_Stats.py` | Dataproc, via WDL | Steps 1 and 2: per-site statistics for one cohort, no filtering |
| `pipeline/Step3_Filter_Rename_Combine.py` | Dataproc, via WDL | Step 3: apply filters, cut both VDSs to kept sites, rename, combine |
| `pipeline/Step4_Densify_Export.py` | Dataproc, via WDL | Step 4: densify and export the VCF (old Step 5b) |
| `pipeline/Eval_Callset.py` | Dataproc, via WDL | unfiltered genotypes of the truth samples |
| `filter_sweep.py` | laptop / Terra terminal | score filter sets offline |
| `check_nonpar_representation.py` | Hail cluster terminal | how each VDS encodes chrX non-PAR |
| `infer_sex_from_vds.py` | Hail cluster terminal | make a sex table for a cohort that has none (AnVIL) |
| `vcfdist_report_sex.py` | laptop / Terra terminal | summarise a VcfdistEvaluation run (patched copy) |
| `export_run_site_lists.py` | Hail cluster terminal | optional: site lists for the report's missed-variant section |

`BeaglePhasing.wdl` and `VcfdistEvaluation.wdl` are used unchanged.

## 1. One-time preparation

1. Upload the four `pipeline/*.py` scripts to the workspace bucket, e.g. `gs://<bucket>/scripts/v2/`.
2. Import `hail_in_wdl_optimized.wdl` into Terra as a new method (or a new snapshot).
3. Sex tables, one per cohort, in the bucket: two columns, sample then sex, tab-, comma- or
   space-separated. `XX`/`F`/`FEMALE` and `XY`/`M`/`MALE` are recognised; every other label
   (`X0`, `XXY`, `NA`, `ambiguous`, ...) is ignored, and a header line is harmless.
   The AnVIL table must use the sample IDs as they are in the AnVIL VDS (no `_C` suffix).

   AnVIL has no sex table yet, so make one from the data, on a Hail cluster terminal:
   ```
   # 3a. check the method on AoU, where the answer is known (prints given label x inferred class)
   python3 infer_sex_from_vds.py \
     --vds gs://prod-drc-broad/v8/wgs/vds/aou_srwgs_short_variants_v8r1.vds/merged \
     --compare-sex-tsv gs://<bucket>/.../aou.foxtrot.sample_sex.tsv \
     --out gs://<bucket>/tmp/aou.inferred_sex.tsv --temp-dir gs://<bucket>/tmp/infer_sex_aou
   # 3b. make the AnVIL table
   python3 infer_sex_from_vds.py \
     --vds gs://fc-secure-759ddb91-68ac-4234-8fab-9c9325c3484f/CCDG_consent.vds --xy-max-het 0.3 \
     --out gs://<bucket>/resources/anvil.inferred_sex.tsv --temp-dir gs://<bucket>/tmp/infer_sex_anvil
   ```
   It reads only variant data, from 100 short windows spread along chrX non-PAR and 20 on chr20,
   and combines two signals per sample: the chrX / autosome depth ratio at non-reference calls
   (about 1 with two X copies, about 0.5 with one) and the share of chrX non-reference calls that
   are heterozygous. A sample is `XX` or `XY` only when the two agree; otherwise `ambiguous`.
   Check the three printed histograms: the depth ratio should have two separate peaks, the
   two-copy group one high het peak, and the one-copy group a peak at 0. In 3a, the cross-table
   should put `XX` in `XX` and `XY`, `X0`, `XYY` in `XY`, with few `ambiguous`. AoU males are
   recognised by depth and by their haploid calls; the low-het rule that AnVIL males depend on is
   exercised in AoU by the diploid-encoded `X0` and `XYY` samples. Use the output of 3b as
   `anvil_sex_tsv_url`. (`--xy-max-het 0.3` for AnVIL: its diploid-encoded males have a het
   fraction of up to about 0.25, against 0 for AoU's haploid-encoded males.)
4. Region BED files (UCSC, 0-based start):
   - test region `chrX_5Mb.bed`: `chrX<tab>0<tab>5000000`
   - smoke test `chrX_200kb.bed`: `chrX<tab>2680000<tab>2880000` (spans the PAR1 boundary)
5. Recommended once: on a Hail cluster terminal, run `check_nonpar_representation.py` for the
   AnVIL VDS with its new sex table (usage in the file header; AoU non-PAR is already done). It
   confirms the sex table matches the sample IDs and shows how often AnVIL males are called het.
   Also run it on AoU for a PAR window (`--interval chrX:1000000-1020000`) to confirm that XY calls
   are diploid there: the pipeline tests all samples in PAR on that assumption.

## 2. Smoke test (about 200 kb; do this first)

Run `hail_in_wdl_optimized.wdl` once with everything switched on, into its own folder:

| Input | Value |
|---|---|
| `aou_vds_url` | `gs://prod-drc-broad/v8/wgs/vds/aou_srwgs_short_variants_v8r1.vds/merged` |
| `anvil_vds_url` | `gs://fc-secure-759ddb91-68ac-4234-8fab-9c9325c3484f/CCDG_consent.vds` |
| `bed_url` | `gs://<bucket>/beds/chrX_200kb.bed` |
| `aou_sex_tsv_url`, `anvil_sex_tsv_url` | the two sex tables |
| `contigs` | `["chrX"]` |
| `prefix` | short lowercase label, e.g. `smoke` |
| `gcs_project` | the workspace's Google project |
| `submission_script_site_stats` | `.../Step1_2_Site_Stats.py` |
| `submission_script_step3` | `.../Step3_Filter_Rename_Combine.py` |
| `submission_script_step4` | `.../Step4_Densify_Export.py` |
| `submission_script_eval_callset` | `.../Eval_Callset.py` |
| `output_bucket_path` | `gs://<bucket>/chrX_smoke_<date>` (a folder not used before) |
| `run_site_stats` | `true` |
| `run_panel` | `true` |
| `export_aou_only` | `true` |
| `eval_sample_names` | the truth sample IDs (the 10 used for vcfdist) |
| `ref_fasta`, `ref_fasta_index` | GRCh38 FASTA and `.fai` (the ones `BeaglePhasing.wdl` uses) |

What to check in the Dataproc job logs (each step prints these):
- Steps 1 and 2: `N samples, N XX, N XY, N without a usable label`. A large unlabelled count means
  the sex table IDs do not match. `variant partitions / reference partitions after interval filter`:
  note these numbers; if they are far below the cluster's core count, the region is being read in
  too few tasks and will not parallelise.
- Step 3: `AoU sites / pass / AnVIL excluded loci / kept` and the resolved thresholds.
- Step 4: the combined VDS partition counts, and `XY labels: N AoU, N AnVIL ...; N of N samples in
  this VDS matched`. The matched count should be close to the two label counts added together; if
  the AnVIL part is missing, the `_C` suffix or the AnVIL sex table is wrong.

Then confirm the outputs exist: the two `site_stats.tsv.bgz`, the two final VCFs with `.tbi`, and
`eval_callset_vcf`. If a step fails, fix it here; the 5 Mb run uses exactly the same code.

### Rerunning only Steps 3 and 4 (after a failed or aborted Step 3 or 4)

1. If the workflow is still running, abort it in Terra. Aborting kills the small VM that would have
   deleted the Dataproc cluster, so delete the cluster yourself (Step 3's lifetime is seven days):
   ```
   gcloud dataproc clusters list   --project <gcs_project> --region us-central1
   gcloud dataproc clusters delete <prefix>-chrx-hail-step3-... --project <gcs_project> --region us-central1
   ```
2. Upload the corrected script(s) over the old ones in the scripts folder.
3. Resubmit with the same inputs except:

   | Input | Value |
   |---|---|
   | `run_site_stats` | `false` (reuses the statistics already under `output_bucket_path`) |
   | `run_panel` | `true` |
   | `eval_sample_names` | `[]` |
   | `filter_tag` | a label not used before in this folder, e.g. `t2` |

   The new `filter_tag` is required: the earlier attempt left partial outputs at its paths, and the
   combiner will not write over a VDS it finds there.
4. On 200 kb, Step 3 should take minutes. Its log should show, in order: the resolved thresholds,
   `AoU sites / pass / AnVIL excluded loci / kept`, the two filtered VDS writes, then the combiner.
   If one task of a stage runs far longer than the rest, send the stage's progress line.

## 3. Phase A: statistics and evaluation callset (5 Mb)

Same inputs as the smoke test except:

| Input | Value |
|---|---|
| `bed_url` | `gs://<bucket>/beds/chrX_5Mb.bed` |
| `prefix` | e.g. `chrx5mb` |
| `output_bucket_path` | `gs://<bucket>/chrX_5Mb_<date>` (call this `<OUT>`; all later phases reuse it) |
| `run_panel` | `false` |
| `export_aou_only` | `false` |

This runs Steps 1 and 2 (in parallel) and the evaluation callset. Outputs:
- `<OUT>/chrX.step1.aou.site_stats.ht` and `.tsv.bgz`
- `<OUT>/chrX.step2.anvil.site_stats.ht` and `.tsv.bgz`
- workflow output `eval_callset_vcf` and `eval_callset_vcf_tbi`

These are computed once. Changing filters never needs them again.

## 4. Phase B: choose filters offline

### B1. Score the unfiltered evaluation callset with vcfdist (once)

Run `VcfdistEvaluation.wdl` with the inputs of your earlier evaluations, changing only:

| Input | Value |
|---|---|
| `eval_vcf` | `eval_callset_vcf` from Phase A |
| `eval_vcf_idx` | `eval_callset_vcf_tbi` from Phase A |
| `do_naively_phase` | `true` (the callset is unphased) |

Keep `eval_sample_names`, `truth_sample_names`, `sample_sex_tsv`, `truth_vcfs`,
`confident_regions_bed_files`, `region`, the reference, `vcfdist_bed_files` (PAR1 then non-PAR)
and `labels_per_stratification` as before. Note the submission ID.

### B2. Look at the male het fraction

```
python3 - <<'PY'
import pandas as pd
for name in ("step1.aou", "step2.anvil"):
    d = pd.read_csv("chrX.%s.site_stats.tsv.bgz" % name, sep="\t", compression="gzip", na_values=["NA"])
    d = d[(d.region == "nonPAR") & (d.n_xy_nonref >= 20)]
    print(name, "non-PAR sites with >= 20 male carriers:", len(d))
    print(d.xy_het_frac.describe(percentiles=[.5, .9, .99, .999]))
PY
```
(copy the two TSVs down first with `gcloud storage cp`).

What the 5 Mb tables showed (2.2 Mb of non-PAR):

- **AoU**: almost 0 everywhere (99% of sites below 0.03). Male calls are haploid; the small
  remainder comes from diploid-encoded single-X samples (`X0`, `XYY`) that the inferred sex table
  labels `XY`.
- **AnVIL** (44,670 sites with at least 20 male carriers): two groups. 73% of sites have a
  fraction of 0.2 or less (median about 0.04); 22% have a fraction above 0.7, most of them above
  0.9; only 5% lie between. A threshold anywhere from 0.3 to 0.5 separates them.
- The high group has a median depth at carriers of 12.3 (the depth of an AnVIL male) against
  22.3 for the rest, so its carriers are almost all male and all called het. That is the pattern
  expected from Y-homologous reads mapping to X.
- About half of the high group currently fails the AnVIL depth filter (cut-off 12) by accident.
  **Relaxing AnVIL depth in non-PAR therefore lets these sites back in unless
  `max_xy_het_frac` is set at the same time.** Test the two together (see B3).

`max_xy_het_frac` only looks at sites with at least `xy_het_min_carriers` male carriers (default
20, as in the snippet above): at rarer sites a single stray male het gives a fraction of 1, and
with `mask_xy_nonpar_hets` on, that call is masked anyway. For AnVIL a threshold near the 0.04
baseline would remove most sites; stay in the 0.3 to 0.5 gap.

### B3. Score filter sets

```
python3 filter_sweep.py \
  --aou-stats   <OUT>/chrX.step1.aou.site_stats.tsv.bgz \
  --anvil-stats <OUT>/chrX.step2.anvil.site_stats.tsv.bgz \
  --gcs-prefix  gs://<bucket>/submissions/<B1 submission folder> \
  --vcfdist-dir vcfdist_evalset --labels PAR1,nonPAR \
  --config aou_flag_in_par:aou.excess_het=flag \
  --config par_relaxed:both.par.min_call_rate=0.8,both.par.min_avsad=10,both.par.min_gq=25 \
  --config par_callrate_only:both.par.min_call_rate=0.8 \
  --config anvil_xy_het_30pct:anvil.nonpar.max_xy_het_frac=0.3 \
  --config anvil_xy_het_50pct:anvil.nonpar.max_xy_het_frac=0.5 \
  --config anvil_nonpar_depth6:anvil.nonpar.min_avsad=6 \
  --config anvil_depth6_xyhet50:anvil.nonpar.min_avsad=6,anvil.nonpar.max_xy_het_frac=0.5 \
  --config anvil_depth6_xyhet30:anvil.nonpar.min_avsad=6,anvil.nonpar.max_xy_het_frac=0.3 \
  --config anvil_nodepth_xyhet50:anvil.nonpar.min_avsad=none,anvil.nonpar.max_xy_het_frac=0.5 \
  --config anvil_depth6_gq20_xyhet50:anvil.nonpar.min_avsad=6,anvil.nonpar.min_gq=20,anvil.nonpar.max_xy_het_frac=0.5 \
  --config no_excess_het:both.excess_het_p=none \
  --config excess_het_with_effect:both.excess_het_f=-0.3 \
  --out filter_sweep.tsv
```
The AnVIL non-PAR configurations come in sets so each change can be read alone and combined:
the male-het filter alone, the depth relaxation alone (expected to readmit the all-het sites), and
the two together.

The excess-het test is Hail's one-sided exact test run per alt allele on diploid calls (XX samples
in non-PAR, everyone elsewhere); a site fails on its smallest per-allele p-value. `excess_het_f`
adds an effect-size condition: the site only fails if that allele's inbreeding coefficient is also
below the value given. The statistics tables carry the per-allele counts, p-values and coefficients
(`eh_*` columns) and the site summary (`p_excess_het`, `excess_het_f`).

Add `--earlier-runs` to also score the filters of the two earlier panels, rebuilt from the same
tables: `pre_handoff_run` (original thresholds, AoU's precomputed ExcessHet flag everywhere, the
all-sample test in AnVIL) and `post_report_run_*` (PAR thresholds 10 / 0.8 / 25, AoU flag in PAR and
the all-sample test in AoU non-PAR, the all-sample test in AnVIL). The post-report run has two rows
because the result depends on whether densified XY reference calls in AoU non-PAR counted as
hom-ref; the row whose `sites_kept` is nearer that run's Step 1 site count is the right one.

Each `--config` is `NAME:SPEC`; `SPEC` is exactly what the WDL takes as `filter_spec`. The spec
grammar is at the top of the WDL. `original` (the defaults) is always included.

How to read the output:
- The offset line should show about 1.000 for one of the two offsets. Below 0.9 the script warns:
  the statistics and the callset are not from the same region or run.
- `truth SNP rows ... not in the AoU site table`: truth the AoU VDS does not contain as a matching
  allele. No filter can recover these.
- `(no filters)` is the ceiling: every site in the table kept.
- For each stratum and sex, one row per filter set: mean per-sample TP / FP / FN, precision,
  recall, F1 for SNPs. Compare each row with `original`.

Limits (also in the script header): SNPs only; genotypes are scored before phasing and with naive
phasing, so use this to rank filter sets and confirm the winner with a full run (Phases C to E).
Sites with AoU max alt AC below 2 have no dense statistics, so `aou.min_ac` cannot be swept below 2
without re-running Step 1 with `--dense_min_ac`.

## 5. Phase C: build the panel with the chosen filters

Run `hail_in_wdl_optimized.wdl` again with the Phase A inputs, changing:

| Input | Value |
|---|---|
| `output_bucket_path` | the same `<OUT>` |
| `run_site_stats` | `false` (reuses the Phase A statistics) |
| `run_panel` | `true` |
| `filter_spec` | the winning spec, or `none` for the defaults |
| `filter_tag` | a short label for this filter set, e.g. `f1`. Use a new one for every filter set |
| `export_aou_only` | `true` if you want the AoU-only comparison arm (section 8) |
| `mask_xy_nonpar_hets` | `true` (default): Step 4 writes het calls of XY samples in chrX non-PAR as missing, so Beagle imputes them. `false` keeps the genotypes as called |
| `eval_sample_names` | `[]` |

Outputs (strings holding bucket paths):
- `final_vcf`: `<OUT>/aou-anvil.final.step5b.f1.chrX.vcf.bgz` (+ `.tbi`)
- `aou_only_vcf`: `<OUT>/aou-only.final.step5b.f1.chrX.vcf.bgz` (+ `.tbi`) when requested
- `<OUT>/chrX.step3.f1.aou.pass.ht`, `...anvil.exclusion.ht`, `...keep.ht`

With `filter_tag` left empty and `filter_spec` `none`, the final VCF has the old pipeline's name.

About `mask_xy_nonpar_hets`: a male has one X in non-PAR, so a het call there cannot be right. AoU
males are haploid-encoded and have none. AnVIL males are diploid-encoded, and a typical one has 5
to 10% of his non-PAR non-reference calls heterozygous; left in, Beagle phases them as real and the
male contributes two different X haplotypes. With masking on, those calls are `./.` in the VCF.
AC / AF / AN / HC in INFO are counted after masking; the mean GQ and AVSAD fields are not changed
by it. Samples labelled anything other than `XY` in the sex tables (including `ambiguous`) are
not masked.

## 6. Phase D: phase with Beagle

Run `BeaglePhasing.wdl` unchanged:

| Input | Value |
|---|---|
| `input_vcf` | `final_vcf` from Phase C |
| `ref_fasta`, `ref_fasta_index` | GRCh38 FASTA and `.fai` |
| `run_par1`, `run_nonpar` | `true` |
| `run_par2` | `false` for the 5 Mb region (`true` for all of chrX) |
| `prefix` | include the filter tag, e.g. `aou-anvil.f1` |

Everything else at its defaults. Outputs needed next: `normalized_vcf`, `normalized_vcf_tbi`.
`beagle_gc_logs` shows how much memory Beagle really used, for right-sizing the machine.

## 7. Phase E: evaluate the panel

1. Run `VcfdistEvaluation.wdl` as before with `eval_vcf` = `normalized_vcf`, `eval_vcf_idx` =
   `normalized_vcf_tbi`, and `do_naively_phase` = `false`.
2. Summarise it:
   ```
   python3 vcfdist_report_sex.py \
     --gcs-prefix gs://<bucket>/submissions/<this submission folder> \
     --labels PAR1,nonPAR --out-dir vcfdist_f1
   ```
   `report.txt` section 3 has precision and recall by region and sex; section 6 has switch and flip
   rates. `evaluation_summary.tsv` in the submission is the table you compared before.
3. Optional, for the missed-variant breakdown (section 7 of the report): on a Hail cluster run
   `export_run_site_lists.py --step1-ht <OUT>/chrX.step3.f1.aou.pass.ht --anvil-exclusion-ht
   <OUT>/chrX.step3.f1.anvil.exclusion.ht --interval chrX:1-5000000 --temp-dir ... --out-dir <folder>`
   into a folder that also holds the older `aou.sites.tsv.bgz` and `anvil.sites.tsv.bgz`, then add
   `--site-tables <folder>` to the command above.

To try another filter set, repeat Phases C to E with a new `filter_tag`.

## 8. What does adding AnVIL buy (phasing)?

With `export_aou_only = true` in Phase C you have two VCFs at exactly the same kept sites:
`final_vcf` (AoU + AnVIL) and `aou_only_vcf` (AoU alone).

1. Run Phase D on each (`prefix` `aou-anvil.f1` and `aou-only.f1`). That is two Beagle runs.
   Masking (`mask_xy_nonpar_hets`) applies to both VCFs but in practice only changes AnVIL
   samples, so the AoU-only arm is the same either way.
2. Run Phase E on each (`--out-dir vcfdist_f1_combined`, `vcfdist_f1_aou_only`).
3. Compare, in section 6 of the two `report.txt` files, the table "Pooled over samples": switches,
   flips and opportunities per region and sex.

Read the comparison with care. The earlier 10-sample evaluation had 70 switch errors in PAR1 and 6
in non-PAR (where only XX samples are informative). Differences smaller than roughly a third in PAR1
are within the noise of counts that size, and non-PAR cannot resolve anything. More truth samples
or a longer region is what makes this question answerable. It also measures the gain for AoU
samples only: no AnVIL sample has truth.

## 9. Moving to all of chrX

- `bed_url`: one line, `chrX<tab>0<tab>156040895`. `run_par2 = true` in `BeaglePhasing.wdl`.
- Cluster lifetime: Steps 1, 2 and 4 are deleted after two days (`time_to_live_minutes`), Step 3
  after seven. Raise these, or the worker counts, if the 5 Mb timings extrapolate past them.
- Nothing here scatters within a contig. If one cluster per step is too slow for 156 Mb, the next
  change is to scatter the Hail steps over region BEDs and join the VCFs before Beagle.
- Send me the `*.report` files (hours per step) from the 5 Mb run before sizing this.

### Autoscaling (optional)

By default every Hail step starts a fixed cluster (2 workers + 8 preemptible). To let Dataproc
add and remove preemptible workers instead, create a policy once and name it in the inputs.

```bash
cat > hail-autoscale.yaml <<'YAML'
workerConfig:
  minInstances: 2
  maxInstances: 2
secondaryWorkerConfig:
  minInstances: 0
  maxInstances: 50
basicAlgorithm:
  cooldownPeriod: 120s
  yarnConfig:
    scaleUpFactor: 0.2
    scaleDownFactor: 0.5
    gracefulDecommissionTimeout: 1200s
YAML
gcloud dataproc autoscaling-policies import hail-autoscale-w2-s50 \
    --source hail-autoscale.yaml --project terra-07e91eb0 --region us-central1
```

Then set `"RunAoUAnvilMergeFilterAndQc.autoscaling_policy": "hail-autoscale-w2-s50"`. With it set:

- every Hail step (1, 2, 3, 4, AoU-only 4, evaluation callset) uses the policy, and each task's
  `num_preemptible_workers` is ignored. `""` (the default) keeps the fixed counts;
- keep `num_workers` at the policy's primary count (2 above), and create the policy in the same
  project and region as the workflow's `gcs_project` and `region`;
- `maxInstances` is the cost ceiling. Check the project's preemptible CPU quota covers it
  (50 x n2d-highmem-8 = 400 vCPUs);
- to use it for one step only, leave the workflow input empty and set that task's own input,
  for example `"RunAoUAnvilMergeFilterAndQc.Step3.autoscaling_policy"`.

Not yet run: try it on the 200 kb region first. Two things could fail at cluster start, both
immediately and with a clear message: the workflow's service account lacking permission to use
autoscaling policies, or this Hail version's `hailctl` rejecting the flag.

## 10. If something fails

| Symptom | Likely cause |
|---|---|
| Step 1/2 exits with "No VDS sample matched an XX (or an XY) label" | sex table IDs do not match the VDS (AnVIL: no `_C` suffix) |
| Step 3 exits with "unknown filter ..." | typo in `filter_spec`; it is checked before the cluster does any work |
| Step 3 fails in the combiner saying the output exists | this `filter_tag` was already used in `<OUT>`; pick a new one |
| Step 3 exits with "No site passes the filters" | a threshold is too strict, or `aou.min_ac` was switched off (sites below `dense_min_ac` then have no dense statistics) |
| One task of a Step 3 stage runs for hours while the rest finished | empty reference padding rows reached the combiner; make sure the Step 3 script is the one with `trim_padding` |
| Executors never start, job hangs at 0 tasks | `executor_memory` does not fit `worker_machine_type` |
| A rerun finishes instantly with stale outputs | call caching returned old path strings; turn caching off for that submission |
| `filter_sweep.py` warns that under 90% of matched SNPs were found | statistics and evaluation callset come from different regions or runs |

## Follow-ups after the 5 Mb evaluation (in priority order)

### 1. Re-score on normalised VCFs

Run `VcfdistEvaluation` on the split, left-aligned VCFs of three panels (original filters, final
AoU-only, final AoU+AnVIL), with the same 50 samples, truth VCFs, BEDs and settings. This settles
whether the lower indel recall of AoU+AnVIL against AoU-only is a scoring artifact of unsplit
multiallelic records. Keep the original run's per-variant outputs: they are needed to attribute each
missing truth SNP to the filter that removed it. The interim (post-report) filters are dropped from
the report.

### 2. Excess-het threshold where it matters

Build the statistics (`run_panel=false`) for a few Mb of non-PAR containing X/Y-homologous sequence
(the X-transposed region, roughly chrX:89-93 Mb; check the coordinates). Look at how many sites
`both.excess_het_f=-0.3` removes there, and at the AnVIL male het fraction.

### 3. Multiallelic sites lost to excess-het

Result on the 5 Mb region: 59 sites are removed by the excess-het test alone, and three of them
carry a well-behaved second allele (chrX:457,623, chrX:499,005, chrX:2,827,053). Per-allele masking
was not built. Repeat the check on the region from item 2:

```bash
python3 check_multiallelic_excess_het.py \
  --aou-stats <OUT>/chrX.step1.aou.site_stats.tsv.bgz \
  --anvil-stats <OUT>/chrX.step2.anvil.site_stats.tsv.bgz \
  --filter-spec both.excess_het_f=-0.3,anvil.nonpar.min_avsad=6
```

The script needs `filter_sweep.py` in the same folder. It lists each site removed by the excess-het
test alone, the cohort whose test fails, and any bystander allele (|F| < 0.1, at least 100 hets, not
failing itself). If bystanders are numerous there, the fix is per-allele masking in Step 4: keep the
site, drop the failing allele, and set calls carrying it to missing.

### 4. PAR2

Representation check in both cohorts and a build of the region.

### 5. Scaling plan for all of chrX

At the 5 Mb rate, Step 4 is about 60 h and Step 1 about a day on 2 + 8 workers. Before the first
full run, choose between scattering the Hail steps by interval and larger or autoscaled clusters
(the local SSD quota limits the second).

### 6. A fairer test of what AnVIL adds

The 50 benchmark samples are AoU participants already in the panel. An imputation test on held-out
samples, or benchmark samples from outside AoU, measures the panel's actual use.

### 7. Lower priority

- The 0.35-point offset between the sweep and vcfdist's summary.
- An independent recount of multiallelic excess-het counts.
- The 439 AnVIL samples without a sex label.
