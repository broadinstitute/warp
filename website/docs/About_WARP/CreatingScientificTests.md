# How to Create a Scientific Test

A guide to the **scientific design** of a WARP pipeline test: deciding what "correct"
means for a pipeline and how to check it in CI. It is pipeline-agnostic, drawn from
building the scANVI (SCVI/SCANVI label-transfer) test, and is meant to be handed to an
agent standing up a test for a different pipeline.

## Scope — read this first

This document does **not** cover the file layout, CI YAML, Dockstore registration, or
truth-bucket mechanics. Those are already written down; do not duplicate them:

- **`AGENTS.md` → "Registering a CI test for a pipeline (Plumbing / Scientific)"** — the
  7-artifact checklist (`Test<Pipeline>.wdl`, `Verify<Pipeline>.wdl`, `test_<pipeline>.yml`,
  `.dockstore.yml` publish, truth seeding, womtool blind spots).
- **`AGENTS.md` → "Test inputs"** — where test JSONs live and how stale inputs break silently.
- **[TestingPipelines.md](TestingPipelines.md)** — the canonical framework reference.
- **[VersionAndReleasePipelines.md](VersionAndReleasePipelines.md)** — version/release rules.

Everything below is the part those docs don't cover: the *judgment* that goes into a
scientific test. Follow the mechanics from the links; make the decisions from here.

## What a Scientific test is (and is not)

Every WARP pipeline test runs the pipeline on Terra and compares its outputs to a stored
truth set. The two test *kinds* answer different questions:

| | Plumbing | Scientific |
|---|---|---|
| Question it answers | "Does the pipeline **run and wire together**?" | "Does the pipeline produce a **scientifically correct** result?" |
| Data | Tiny / synthetic / truncated | Realistic, representative of a real use case |
| When | Every PR | PRs to master, or on demand |
| Cost | Seconds–minutes, cheap | Minutes–hours, real compute |
| What a failure means | Something is broken structurally | The science regressed |

The defining test for whether you have a *Scientific* test: **would it catch a result that
runs cleanly but is biologically/statistically wrong?** If shrinking the data or capping
the work would let a wrong-but-runnable result slip through, that shrunk version is a
Plumbing test, not a Scientific one. (Deriving the Plumbing test from the Scientific one is
covered at the end.)

## Step 1 — Define "correct" as a set of invariants

Before touching data or WDL, write down — in plain language — what must be true of a
correct output. Sort each statement into one of two buckets, because they are verified
differently:

- **Structural invariants** — must hold *exactly*, every run. Shape, presence of expected
  fields/columns/files, row counts, no forbidden values. These are cheap and non-negotiable.
- **Scientific/distributional invariants** — hold *approximately*, within a tolerance,
  because the pipeline is stochastic or floating-point sensitive. "The answer is close to
  the reference," not "the answer is byte-identical."

Example (scANVI, a stochastic label-transfer pipeline):

> Structural: output h5ad has the **same number of cells** as truth; the annotation column
> is **present**; the predicted-label vocabulary is a **subset** of truth's (the model may
> not invent labels the reference never had).
> Distributional: **per-cell-type proportions correlate** with truth at or above a
> threshold.

Doing this first turns "compare the outputs" into a concrete, reviewable specification and
tells you exactly what the verification WDL must check.

## Step 2 — Choose test data that can actually expose a regression

Scientific test data should be **realistic and representative, but bounded**:

- **Representative of a real use case.** Use a real input + real reference, at the scale and
  configuration a user would actually run (scANVI: a real reference atlas and a real query
  set at a meaningful taxonomy level). A toy input can't validate science.
- **Bounded for CI.** Big enough to exercise the science, small enough to finish in a CI
  time/cost budget. Subsample rather than fabricate, and subsample in a way that preserves
  the property you're testing (e.g. keep enough cells per class that per-class proportions
  are meaningful).
- **Able to expose the failure you care about.** Ask "what would a wrong result look like?"
  and confirm the data + invariants would catch it. If every plausible bug still passes,
  the test proves nothing.
- **Watch for leakage / trivialization.** For ML pipelines, make sure the test isn't
  trivially memorized (query identical to a training row) unless that's deliberately the
  thing under test.

Host the data under the public test bucket per the AGENTS.md mechanics. Record where it came
from and how it was subsampled (a `test_data_overview.md` next to the inputs, as several
pipelines do), so the next person knows what the numbers *should* look like.

## Step 3 — Design the verification (the hard part)

The verification WDL (`Verify<Pipeline>.wdl`) is where the science lives. Match the
comparison to the pipeline's determinism:

### Deterministic pipelines
If the same inputs always produce the same bytes, compare closely — checksums, exact array
equality, or a tight numeric tolerance. Any drift is a real change.

### Stochastic / ML / float-sensitive pipelines
If outputs vary run-to-run (random seeds, GPU nondeterminism, threading), **exact comparison
is wrong** — it will be flaky and will force people to re-bless truth constantly, destroying
its value. Instead verify the invariants from Step 1:

1. **Check structural invariants exactly.** Cell/row counts, presence of columns/files,
   label vocabulary containment. These are deterministic even when values aren't.
2. **Check distributional invariants against a threshold.** Correlate distributions, compare
   summary statistics, bound a divergence — whatever captures "close to the reference." In
   scANVI this is a correlation of per-cell-type proportions with a `min_proportion_corr`
   threshold.
3. **Handle degenerate cases explicitly.** Single-class inputs, zero-variance vectors, empty
   optional outputs — decide what these *mean* (pass or fail) rather than letting a `NaN`
   correlation silently pass or crash. (scANVI treats an all-identical proportion vector as a
   perfect match, and a `NaN` correlation otherwise as `0.0` so the threshold check fails
   rather than passes on `NaN`.)

### Calibrating the tolerance
A threshold that's too tight is flaky; too loose passes garbage. Calibrate empirically:
generate truth, then run the pipeline **a few more times** and measure the natural
run-to-run variance of your distributional metric. Set the threshold comfortably *below* the
worst honest run but *above* what a real regression would produce. Document the chosen number
and why (scANVI: `min_proportion_corr = 0.95`).

### Fail loudly on things "tolerant" must NOT excuse
Tolerance is about *values*, not *existence*. A whole output disappearing, a column
vanishing, or an output appearing on only one of {test, truth} is a regression even in a
stochastic pipeline — assert on it. scANVI's ATAC output is optional (multiome only), so its
verify **errors** if exactly one of test/truth has it, rather than silently skipping the
comparison.

### Exclude genuinely non-comparable artifacts — but archive them
Some outputs can't be meaningfully compared (e.g. trained model weights, which differ every
run). Don't verify them, but still **copy them to the results/truth bucket** so they're
available for debugging and reuse. Say so in a comment (scANVI archives `scanvi_model_out`
without verifying it).

### Keep the compare task pipeline-local
Put the pipeline's compare task in its own `Verify<Pipeline>.wdl`, **not** the shared
`verification/VerifyTasks.wdl` — editing the shared file triggers every pipeline's CI (see
AGENTS.md). Keep the verifier's own container light and its dependencies minimal; it runs on
every compare.

## Step 4 — Truth is a scientific judgment, not just a file

Seeding truth (`updateTruth: true`) **declares the current outputs correct**. Treat that as a
sign-off, not a mechanical step:

- **Generate truth from a known-good pipeline version** — ideally the last released one, or a
  version whose output you've inspected.
- **Inspect before you bless.** Look at the update-truth run's outputs for scientific sanity
  (do the predicted labels/distributions look right?), not merely that files exist. Truth is
  the yardstick for every future run; a wrong truth silently passes wrong results forever.
- **Re-bless deliberately when outputs legitimately change.** A change that intentionally
  alters outputs (new metadata field, model change, a genuinely better result) *requires*
  updating truth — and requires re-inspecting it. Note the reason in the changelog.
- **Truth is keyed by the test JSON filename**, not `input_id` (see AGENTS.md step 6): a new
  test case needs its truth seeded first, or the first compare run fails with nothing to diff
  — the most common "new test" failure.

## Step 5 — Bound the cost honestly

Scientific tests use real compute, so right-size but don't cripple:

- Use a modest but sufficient machine (scANVI verifies on a single T4, not a large multi-GPU
  box).
- Do **not** cap the scientific work in a way that changes the science (e.g. slashing
  training epochs) — that's what the Plumbing test is for. If you find yourself wanting to,
  you're building a Plumbing test.
- Lean on `useCallCache` during development to avoid re-running unchanged upstream steps.

## The development loop

Standing up a Scientific test is iterative:

1. Write inputs + `Test<Pipeline>.wdl` + `Verify<Pipeline>.wdl` (mechanics per AGENTS.md);
   validate with womtool.
2. **Seed truth** (`updateTruth: true`) and **inspect** the outputs (Step 4).
3. Run a **compare** run (`updateTruth: false`). Expect to tune: threshold too tight → tighten
   data or loosen threshold with justification; passes obviously-wrong output → strengthen an
   invariant.
4. Repeat 2–3 until the test reliably passes on good code and you're convinced it would fail
   on a real regression.

## Checklist

- [ ] Invariants written down, split into structural (exact) vs distributional (tolerant).
- [ ] Test data is realistic, bounded, and can expose a real regression; provenance documented.
- [ ] Verification checks structural invariants exactly and distributional ones against a
      *calibrated, documented* threshold.
- [ ] Degenerate cases (single class, NaN, empty optional outputs) handled explicitly.
- [ ] Presence/absence and shape regressions fail loudly (not excused by tolerance).
- [ ] Non-comparable artifacts excluded from comparison but still archived.
- [ ] Compare task lives in `Verify<Pipeline>.wdl`, not the shared `VerifyTasks.wdl`.
- [ ] Truth seeded from a known-good version and visually inspected before blessing.
- [ ] Mechanics (CI YAML, Dockstore publish, permissions) done per AGENTS.md.

---

## Appendix: Creating a Plumbing test from a Scientific test

A Plumbing test is a **shrunk Scientific test that proves the pipeline wires together
end-to-end**, cheaply, on every PR. Reuse the same `Test<Pipeline>.wdl` and
`Verify<Pipeline>.wdl` — only the **input JSON** and a few knobs change:

- **Shrink the data** to the smallest input that still exercises every task (tiny/truncated
  reference and query).
- **Turn cost knobs down** using inputs the pipeline already exposes: cap iterations (scANVI
  sets `max_epochs: 2`), drop the GPU where possible (scANVI's pretrained-model case runs
  `gpu_count: 0` on the CPU task), and shrink mem/cpu/disk.
- **Keep the same verification, but expect looser signal.** With trivial data the
  distributional numbers are meaningless, so the value is mostly the *structural* invariants
  (it ran, outputs exist with the right shape/columns). The tolerant compare still runs; it
  just isn't proving much scientifically — which is fine, that's the Scientific test's job.
- **Select the kind at dispatch:** the CI wrapper picks Plumbing vs Scientific via `testType`
  (and the branch-derived default). Each kind has its **own truth**, seeded separately —
  seed the Plumbing truth too.

Rule of thumb: if a change to the input JSON is all it takes, you're deriving a Plumbing test
correctly. If you need to weaken the *verification* to make trivial data pass, prefer relying
on the structural invariants instead — don't lower the scientific thresholds that protect the
Scientific test.
