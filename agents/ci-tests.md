# CI tests: agent guide

Read this when adding or changing a Plumbing/Scientific test, test inputs, a `Verify*` compare task, or a `.github/workflows/test_*.yml`. The general rules are in [AGENTS.md](../AGENTS.md); the canonical framework reference is [TestingPipelines.md](../website/docs/About_WARP/TestingPipelines.md), and the judgment that goes into a Scientific test is in [CreatingScientificTests.md](../website/docs/About_WARP/CreatingScientificTests.md).

## Test inputs

Test JSON inputs live at two locations per pipeline:

- `pipelines/wdl/<name>/test_inputs/{Plumbing,Scientific}/` — CI test inputs (Plumbing = fast/cheap, Scientific = end-to-end accuracy)
- `pipelines/wdl/<name>/example_inputs/` — example inputs for documentation

When adding a new test input for a regression case, also wire it into the relevant GitHub Actions workflow under `.github/workflows/` so CI picks it up.

**Stale inputs break silently** — they are not validated by womtool. When **adding**, renaming, or removing/renaming inputs, audit `example_inputs/*.json` and `test_inputs/**/*.json` **and** thread the change through the `Test<Name>.wdl` wrapper (see [Registering a CI test](#registering-a-ci-test-for-a-pipeline-plumbing--scientific) step 2 — *forward every testable input*). Mind the framework's `UpdateTestInputs.py` key rewrite: a test-JSON key with more than two dotted parts (e.g. `Pipeline.Sub.input`, reaching into a sub-workflow) is rewritten to `TestPipeline.Pipeline.Sub.input`, which the wrapper rejects as an extra input. Prefer exposing the input at the pipeline's **top level** (`Pipeline.input`) and setting it that way — a two-part key rewrites cleanly to `TestPipeline.input`.

## Registering a CI test for a pipeline (Plumbing / Scientific)

A WARP "Plumbing" or "Scientific" test is a **CI-integrated test that runs the pipeline on Terra and compares its outputs to a stored truth set** — not a local demonstration or a notebook. Standing one up for pipeline `<Pipeline>` requires this whole set of artifacts; a missing piece typically fails only at CI run time (womtool does not catch most of it). Canonical reference: [TestingPipelines.md](../website/docs/About_WARP/TestingPipelines.md). There is also an internal SOP, *"How to Port a Pipeline Test to the New Testing Framework"* (Optimus example) — useful for the Dockstore-publish click-path, but its file paths predate the directory flatten (it shows `pipelines/skylab/`, `tasks/broad/`; use today's `pipelines/wdl/`, `tasks/wdl/`), and its "TestX.wdl modifications" section is Vault→Terra migration only, not relevant to a net-new test. The fastest way to author the Actions YAML and the test wrapper is to copy an existing pipeline's `test_<name>.yml` / `Test<Name>.wdl` and adapt.

1. **Test input JSON(s)** — `pipelines/wdl/<name>/test_inputs/{Plumbing,Scientific}/<Case>.json`, keyed by the **pipeline's** input names (`<Pipeline>.foo`). Host the input data under `gs://pd-test-storage-public/<Pipeline>/input/{plumbing,scientific}/...`. Plumbing = tiny/fast/cheap (every PR); Scientific = realistic end-to-end (PRs to master, or on demand).
2. **Test wrapper** — `verification/test-wdls/Test<Pipeline>.wdl`: imports the pipeline, `Verify<Pipeline>`, `Utilities`, and `TerraCopyFilesFromCloudToCloud`; declares **every pipeline input a test JSON might set** plus the framework-injected `truth_path`, `results_path`, `update_truth` (and `cloud_provider` if the pipeline takes it); sets `meta { allowNestedInputs: true }`. It calls the pipeline, gathers outputs into an `Array[String]`, copies them to `results_path`, copies to `truth_path` when `update_truth`, else runs `GetValidationInputs` and calls `Verify<Pipeline>`. **Forward every testable input** — an input the pipeline has but the wrapper omits crashes the moment a test JSON sets it (see [Test inputs](#test-inputs) above).
3. **Verification** — `verification/Verify<Pipeline>.wdl` compares each output to truth (tolerantly). Keep pipeline-specific compare tasks **here, not in the shared `verification/VerifyTasks.wdl`**, so editing them doesn't trip every other pipeline's path filter. When a `Verify*` comparison fails on tiny float or "1-of-N" differences, it is usually known run-to-run nondeterminism, not a regression — if the randomness comes from a seedable tool, prefer an optional seed input on the pipeline (unset = production behavior) pinned in `Test<Pipeline>.wdl`, then compare exactly (precedent: Optimus `emptydrops_seed`); otherwise compare with a tolerance by extending an existing tolerant task (e.g. `CompareGeneMetricsWithTolerance`, or the tolerant `CompareSAMs` in `VerifyRNAWithUMIs.wdl`) rather than adding a strict check. Any metric absent from a percentage-threshold table (e.g. `CompareAtacLibraryMetrics` in `VerifyTasks.wdl`) defaults to **exact match** — a newly-nondeterministic metric must be added explicitly or it fails on the first sub-unit drift. When the threshold is a **minimum** (e.g. `CompareScanviH5ad`'s proportion-correlation floor), set it below the worst observed across **all** variants — the GEX and ATAC label paths drift differently (observed ~0.875 GEX vs ~0.783 ATAC), so a floor tuned to one flakes the other.
4. **GitHub Actions entry point ("the buttons")** — `.github/workflows/test_<pipeline>.yml`: `on: pull_request` with a `paths:` filter scoped to the pipeline's files (pipeline dir, its tasks, `Verify<Pipeline>.wdl`, `Test<Pipeline>.wdl`, `TerraCopyFilesFromCloudToCloud.wdl`, the two workflow files, `firecloud_api.py`) **and** `workflow_dispatch` with inputs `useCallCache`, `updateTruth`, `testType` (choice Plumbing/Scientific), `truthBranch`. The job `uses: ./.github/workflows/warp_test_workflow.yml` with `pipeline_name: Test<Pipeline>`, `dockstore_pipeline_name: <Pipeline>`, `pipeline_dir`, those inputs, and `secrets: { PDT_TESTER_SA_B64, DOCKSTORE_TOKEN }`; `permissions: { contents: read, id-token: write, actions: write }`. (`warp_test_workflow.yml` is the shared reusable workflow: it creates the Terra method config, submits, polls, copies results, and runs verification when not updating truth.)
5. **Dockstore registration *and publish*** — add a `Test<Pipeline>` entry (`name`, `subclass: WDL`, `primaryDescriptorPath: /verification/test-wdls/Test<Pipeline>.wdl`) to `.dockstore.yml` (the pipeline itself is usually already registered). Then **publish it in the Dockstore UI** — a manual step that's easy to forget and not caught by anything in-repo: wait ~5 min for Dockstore to ingest the new descriptor, go to the [Dockstore dashboard](https://dockstore.org/dashboard), find `warp/Test<Pipeline>` via the *Search Workflows* field, then **Versions → Actions → Set as Default Version → Publish**. The CI resolves the workflow by `dockstore_pipeline_name`, so an unpublished / no-default-version workflow makes the Terra method-config step fail.
6. **Seed truth FIRST** — a brand-new test has no golden files. Run the Actions workflow manually (`workflow_dispatch`) with `updateTruth: true` and the right `truthBranch` to populate the truth bucket; only then do compare-runs pass. A compare-run before truth exists fails with nothing to diff. The truth path is keyed by the **test JSON's filename**, not `input_id`: `gs://pd-test-storage-public/<Pipeline>/truth/<plumbing|scientific>/<truthBranch>/<JSON-basename>/` (e.g. `mouse_v4_snRNA_example.json` → `.../mouse_v4_snRNA_example/`). So changing only `input_id` overwrites the *same* truth key; renaming the JSON file starts a fresh one.
7. **Validate** every new/modified WDL with womtool. Note what womtool does **not** check: the Actions YAML, the JSON keys, the `.dockstore.yml` entry, and whether truth exists — all of which only surface at CI run time.

If a pipeline ships only one test kind, pin it in the caller (e.g. scANVI is Scientific-only → `test_type: ${{ github.event.inputs.testType || 'Scientific' }}`, and `warp_test_workflow.yml` honors an explicit caller override over the branch-derived default).

Test generation utilities live in `verification/test-wdls/scripts/`.

## Sizing Optimus test data

Shrink Plumbing/Scientific inputs by subsetting reads to a genomic **region** (e.g. one chromosome), which keeps the called-cell count high. Do **not** subsample down to a handful of **barcodes**: the library-level doublet step sets its KNN neighbor count to `k = round(min(100, n_called_cells × 0.01))`, so below ~50 cells `k` becomes `0` and the run dies with sklearn's `n_neighbors ... Got 0` error. Keep ≥~500 called cells. See the [Minimum cell count warning](../website/docs/Pipelines/Optimus_Pipeline/README.md).

## Tightening a shared Verify task

Before making a shared compare task stricter (e.g. `CompareH5adFilesGEX` in `verification/VerifyTasks.wdl`, called by six `Verify*.wdl` files), check every caller's `Test<Pipeline>.wdl` for inputs that can route a test onto the tightened path, not just the current test JSONs. Example: Multiome/PairedTag/SlideTags tests all run `counting_mode = "sn_rna"`, but their wrappers accept `"sc_rna"`, which runs EmptyDrops.

## CI / GitHub workflows

The branch model (develop/staging/master) and how the smart test system selects which pipelines to test are documented in [TestingPipelines.md](../website/docs/About_WARP/TestingPipelines.md). Operational specifics:

- **Path-based triggers** — GitHub Actions workflows trigger on the paths a pipeline depends on:

  ```yaml
  paths:
    - 'pipelines/wdl/optimus/**'
    - 'tasks/wdl/StarAlign.wdl'
    - 'tasks/wdl/Metrics.wdl'
    - 'tasks/wdl/Utilities.wdl'
  ```

  When adding a regression test input, wire its path into the relevant workflow so CI picks it up.
- **Cancel scope** — `warp_test_workflow.yml` runs in concurrency group `warp-test-<pipeline_name>-<branch>` with `cancel-in-progress`, so a push to a PR branch or a manual dispatch of the same pipeline on that branch cancels that pipeline's running test there (e.g. a long manual Scientific or truth run). Re-running *other* pipelines' failed jobs is safe.
- **Release builds** — `./scripts/build_pipeline_release.sh -w pipeline.wdl -v version -o output_dir`.

## Stale CI path filters

Each `.github/workflows/test_*.yml` has a `paths:` filter enumerating the WDLs whose changes should trigger that pipeline's test. The filter is maintained by hand, so it drifts out of sync with what the Test WDL actually imports — a watched file gets renamed/deleted (the filter then watches a dead or *wrong* file), or a new sub-workflow/task import is added but never watched. Either way the test silently stops firing on the edits that matter. Scan for both:

```bash
scripts/find_stale_ci_path_filters.sh
```

It reports `DEAD` (watched path missing / points at the wrong file) and `BLIND` (imported by the test but not covered by any watched path or `dir/**` glob) per workflow, and flags any workflow that doesn't resolve to exactly one Test WDL.

**These are CI-only changes** — editing a workflow's `paths:` filter touches no WDL, so there's no version bump or changelog cascade. Fix them anytime; they don't need to ride a pipeline's next real change the way [unused imports](wdl-authoring.md#unused-imports) do.

**Future CI wiring (not done — deliberately):** same rationale as the [unused-import scan](wdl-authoring.md#unused-imports) — the test system is flaky, so this stays a manual scan. When CI is trusted, wire it as an **advisory** (`continue-on-error`) job before ever making it blocking.
