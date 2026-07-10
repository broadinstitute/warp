# Agent Instructions

Single entry point for all agent-driven and AI-coding work in this repository. Start here. This file holds the agent-operational guidance; for the *rules* and *rationale* behind it, follow the links to the canonical human-facing docs below — those are the source of truth, so do not restate them here.

## Authoritative References

| Topic | Document |
| --- | --- |
| Changelog format and language (entry structure, active voice, sample entries) | [website/docs/contribution/contribute_to_warp/changelog_style.md](website/docs/contribution/contribute_to_warp/changelog_style.md) |
| Semantic-versioning classification (major/minor/patch) and release packaging | [website/docs/About_WARP/VersionAndReleasePipelines.md](website/docs/About_WARP/VersionAndReleasePipelines.md) |
| Pipeline best practices (automation, testability, portability, scaling, maintainability) | [website/docs/About_WARP/BestPractices.md](website/docs/About_WARP/BestPractices.md) |
| Requirements a released pipeline must meet | [website/docs/About_WARP/PipelineRequirements.md](website/docs/About_WARP/PipelineRequirements.md) |
| Testing: branch model (develop/staging/master) and how the smart test system selects pipelines | [website/docs/About_WARP/TestingPipelines.md](website/docs/About_WARP/TestingPipelines.md) |
| Docusaurus markdown style (admonitions, tables, tabs, code-block highlighting, cross-refs) | [website/docs/contribution/contribute_to_warp_docs/doc_style.md](website/docs/contribution/contribute_to_warp_docs/doc_style.md) |
| How to build / serve / validate the docs site | [website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md](website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md) |

The human-facing docs are canonical for rules and rationale. Anything *agent-specific* or *operational* lives in this file — do not duplicate human-facing content here; link to it instead.

## Quick Checklist Before Completing Any WARP Task

1. **Validate all modified WDLs** — and every WDL that imports them (transitively). See [WDL Validation (MANDATORY)](#wdl-validation-mandatory).
2. **Changelog and versioning** — follow the [agent automation rules](#changelog-and-versioning) below; for entry format see the [changelog style guide](website/docs/contribution/contribute_to_warp/changelog_style.md), and for what makes a change major/minor/patch see [VersionAndReleasePipelines.md](website/docs/About_WARP/VersionAndReleasePipelines.md).
3. **Cascading version bumps** — when modifying a shared task (`tasks/wdl/`), every pipeline that imports it (directly or transitively) needs either a patch bump + changelog note **or** an explicit "no functional impact" entry. See [Cascading version bumps](#cascading-version-bumps). Validation passing is not sufficient — changelogs must also be updated.
4. **After merging develop**, grep for duplicate `pipeline_version` lines — keep the higher one. See [Git merge conflict resolution](#git-merge-conflict-resolution).
5. **Sub-workflow contract** — removing an input from a shared WDL requires removing it from every caller. See [Sub-workflow input contract](#sub-workflow-input-contract).
6. **Stale example/test inputs** — when you rename a workflow or remove/rename inputs, audit `pipelines/wdl/<name>/example_inputs/*.json` and `test_inputs/**/*.json`; they break silently because they are not checked by womtool.
7. **Touching a pipeline's interface** — also update the pipeline's docs page under `website/docs/Pipelines/<Name>_Pipeline/README.md` and run `yarn --cwd=website build` to catch broken links.

## Repository Structure

WARP is a collection of cloud-optimized WDL (Workflow Description Language) pipelines for biological data processing. It uses a **flattened directory structure**:

- **`pipelines/wdl/`** — all workflow definitions, organized by category (`arrays/`, `atac/`, `dna_seq/`, `genotyping/`, `multiome/`, `optimus/`, `peak_calling/`, `reprocessing/`, `rna_seq/`, …)
- **`tasks/wdl/`** — all reusable WDL tasks (e.g. `Alignment.wdl`, `BamProcessing.wdl`, `FastqProcessing.wdl`, `StarAlign.wdl`, `Utilities.wdl`)
- **`structs/`** — WDL struct definitions for type safety
- **`verification/`** — test workflows that validate pipeline outputs (`verification/test-wdls/` for test implementations)
- **`scripts/`** — build and validation automation

> **History:** Content was previously split between `pipelines/broad/` (DNA-seq, arrays, reprocessing) and `pipelines/skylab/` (single-cell RNA-seq, ATAC-seq, multiome)--any remnants you should offer to clean up. Both are now unified under `pipelines/wdl/` and `tasks/wdl/`. Docker images are maintained separately in [warp-tools](https://github.com/broadinstitute/warp-tools).

## WDL Development and Validation

### Import structure

All WDL files use relative imports:

```wdl
import "../../../tasks/wdl/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../structs/dna_seq/DNASeqStructs.wdl"
```

Import path errors are the most common validation failure.

### Sub-workflow input contract

A WDL file defines a **single input contract for all callers**. You cannot expose an input to one workflow but hide it from another that imports the same WDL. To remove a parameter from one consumer, you must remove it from the shared task/sub-workflow **and** from every caller.

### WDL Validation (MANDATORY)

Womtool validation is REQUIRED for all WDL changes, and must pass for every modified file **and its dependents**. The womtool JAR is expected at a user-specific location such as `~/womtool-90.jar` (with a compatible Java available); a convenience wrapper exists at `scripts/validate_wdls.sh`.

Validate a single file:

```bash
java -jar ~/womtool-90.jar validate path/to/pipeline.wdl
```

**After modifying a task, validate every workflow that imports it** — a task can validate in isolation but break a caller's argument list. Use a loop:

```bash
for f in tasks/wdl/StarAlign.wdl tasks/wdl/CheckInputs.wdl \
          pipelines/wdl/optimus/Optimus.wdl \
          pipelines/wdl/multiome/Multiome.wdl \
          pipelines/wdl/paired_tag/PairedTag.wdl \
          pipelines/wdl/slidetags/SlideTags.wdl \
          pipelines/wdl/slideseq/SlideSeq.wdl \
          pipelines/wdl/smartseq2_single_nucleus_multisample/MultiSampleSmartSeq2SingleNucleus.wdl \
          verification/test-wdls/TestOptimus.wdl \
          verification/test-wdls/TestMultiome.wdl \
          verification/test-wdls/TestPairedTag.wdl \
          verification/test-wdls/TestSlideTags.wdl \
          verification/test-wdls/TestSlideSeq.wdl; do
  printf '%-90s ' "$f"
  java -jar ~/womtool-90.jar validate "$f" 2>&1 | tail -1
done
```

### Test inputs

Test JSON inputs live at two locations per pipeline:

- `pipelines/wdl/<name>/test_inputs/{Plumbing,Scientific}/` — CI test inputs (Plumbing = fast/cheap, Scientific = end-to-end accuracy)
- `pipelines/wdl/<name>/example_inputs/` — example inputs for documentation

When adding a new test input for a regression case, also wire it into the relevant GitHub Actions workflow under `.github/workflows/` so CI picks it up.

**Stale inputs break silently** — they are not validated by womtool. When renaming a workflow or removing/renaming inputs, audit `example_inputs/*.json` and `test_inputs/**/*.json`.

### Registering a CI test for a pipeline (Plumbing / Scientific)

A WARP "Plumbing" or "Scientific" test is a **CI-integrated test that runs the pipeline on Terra and compares its outputs to a stored truth set** — not a local demonstration or a notebook. Standing one up for pipeline `<Pipeline>` requires this whole set of artifacts; a missing piece typically fails only at CI run time (womtool does not catch most of it). Canonical reference: [TestingPipelines.md](website/docs/About_WARP/TestingPipelines.md). There is also an internal SOP, *"How to Port a Pipeline Test to the New Testing Framework"* (Optimus example) — useful for the Dockstore-publish click-path, but its file paths predate the directory flatten (it shows `pipelines/skylab/`, `tasks/broad/`; use today's `pipelines/wdl/`, `tasks/wdl/`), and its "TestX.wdl modifications" section is Vault→Terra migration only, not relevant to a net-new test. The fastest way to author the Actions YAML and the test wrapper is to copy an existing pipeline's `test_<name>.yml` / `Test<Name>.wdl` and adapt.

1. **Test input JSON(s)** — `pipelines/wdl/<name>/test_inputs/{Plumbing,Scientific}/<Case>.json`, keyed by the **pipeline's** input names (`<Pipeline>.foo`). Host the input data under `gs://pd-test-storage-public/<Pipeline>/input/{plumbing,scientific}/...`. Plumbing = tiny/fast/cheap (every PR); Scientific = realistic end-to-end (PRs to master, or on demand).
2. **Test wrapper** — `verification/test-wdls/Test<Pipeline>.wdl`: imports the pipeline, `Verify<Pipeline>`, `Utilities`, and `TerraCopyFilesFromCloudToCloud`; declares **every pipeline input a test JSON might set** plus the framework-injected `truth_path`, `results_path`, `update_truth` (and `cloud_provider` if the pipeline takes it); sets `meta { allowNestedInputs: true }`. It calls the pipeline, gathers outputs into an `Array[String]`, copies them to `results_path`, copies to `truth_path` when `update_truth`, else runs `GetValidationInputs` and calls `Verify<Pipeline>`. **Forward every testable input** — an input the pipeline has but the wrapper omits crashes the moment a test JSON sets it (see *Test inputs* above).
3. **Verification** — `verification/Verify<Pipeline>.wdl` compares each output to truth (tolerantly). Keep pipeline-specific compare tasks **here, not in the shared `verification/VerifyTasks.wdl`**, so editing them doesn't trip every other pipeline's path filter.
4. **GitHub Actions entry point ("the buttons")** — `.github/workflows/test_<pipeline>.yml`: `on: pull_request` with a `paths:` filter scoped to the pipeline's files (pipeline dir, its tasks, `Verify<Pipeline>.wdl`, `Test<Pipeline>.wdl`, `TerraCopyFilesFromCloudToCloud.wdl`, the two workflow files, `firecloud_api.py`) **and** `workflow_dispatch` with inputs `useCallCache`, `updateTruth`, `testType` (choice Plumbing/Scientific), `truthBranch`. The job `uses: ./.github/workflows/warp_test_workflow.yml` with `pipeline_name: Test<Pipeline>`, `dockstore_pipeline_name: <Pipeline>`, `pipeline_dir`, those inputs, and `secrets: { PDT_TESTER_SA_B64, DOCKSTORE_TOKEN }`; `permissions: { contents: read, id-token: write, actions: write }`. (`warp_test_workflow.yml` is the shared reusable workflow: it creates the Terra method config, submits, polls, copies results, and runs verification when not updating truth.)
5. **Dockstore registration *and publish*** — add a `Test<Pipeline>` entry (`name`, `subclass: WDL`, `primaryDescriptorPath: /verification/test-wdls/Test<Pipeline>.wdl`) to `.dockstore.yml` (the pipeline itself is usually already registered). Then **publish it in the Dockstore UI** — a manual step that's easy to forget and not caught by anything in-repo: wait ~5 min for Dockstore to ingest the new descriptor, go to the [Dockstore dashboard](https://dockstore.org/dashboard), find `warp/Test<Pipeline>` via the *Search Workflows* field, then **Versions → Actions → Set as Default Version → Publish**. The CI resolves the workflow by `dockstore_pipeline_name`, so an unpublished / no-default-version workflow makes the Terra method-config step fail.
6. **Seed truth FIRST** — a brand-new test has no golden files. Run the Actions workflow manually (`workflow_dispatch`) with `updateTruth: true` and the right `truthBranch` to populate the truth bucket; only then do compare-runs pass. A compare-run before truth exists fails with nothing to diff. The truth path is keyed by the **test JSON's filename**, not `input_id`: `gs://pd-test-storage-public/<Pipeline>/truth/<plumbing|scientific>/<truthBranch>/<JSON-basename>/` (e.g. `mouse_v4_snRNA_example.json` → `.../mouse_v4_snRNA_example/`). So changing only `input_id` overwrites the *same* truth key; renaming the JSON file starts a fresh one.
7. **Validate** every new/modified WDL with womtool. Note what womtool does **not** check: the Actions YAML, the JSON keys, the `.dockstore.yml` entry, and whether truth exists — all of which only surface at CI run time.

If a pipeline ships only one test kind, pin it in the caller (e.g. scANVI is Scientific-only → `test_type: ${{ github.event.inputs.testType || 'Scientific' }}`, and `warp_test_workflow.yml` honors an explicit caller override over the branch-derived default).

## Changelog and Versioning

Every pipeline carries a `String pipeline_version = "major.minor.patch"` and a cumulative `<PipelineName>.changelog.md`. For entry **format and language** follow the [changelog style guide](website/docs/contribution/contribute_to_warp/changelog_style.md); for what makes a change **major / minor / patch** follow [VersionAndReleasePipelines.md](website/docs/About_WARP/VersionAndReleasePipelines.md). The agent-operational rules below are not in those docs:

1. **Version bumps are once per branch/PR.** If the top changelog entry on the current branch already has an unreleased version bump (i.e. its version is higher than what is in `develop`), append new bullet points to that entry — do not create a new version section.
2. **`pipeline_versions.txt`** must be updated only when a pipeline version number actually changes — not on every commit. Keep it in sync with the WDL `pipeline_version` string and the top of the changelog.
3. **Linked versions** — some pipelines reference a sub-workflow's version (e.g. ImputationBeagle references the ArrayImputationQuotaConsumed version); update both together.

### Cascading version bumps

When a shared task changes, every workflow that imports it may need a version bump. Either bump the version with a changelog note **or** add a bullet explicitly stating there is no functional impact on that pipeline. The concrete dependency chains:

- `tasks/wdl/StarAlign.wdl` → **Optimus**, **SlideSeq**, **Multiome** (via Optimus), **PairedTag** (via Optimus), **SlideTags** (via Optimus), **MultiSampleSmartSeq2SingleNucleus**
- `tasks/wdl/CheckInputs.wdl` → **Optimus** (via `checkOptimusInput`) and its wrappers (**Multiome**, **PairedTag**, **SlideTags**); also **MultiSampleSmartSeq2SingleNucleus** (via `checkInputArrays`, a separate task in the same file). Note: SlideSeq imports this file but never calls any of its tasks — changes to CheckInputs.wdl task signatures do not functionally affect SlideSeq through this import.
- `tasks/wdl/Metrics.wdl` → Most Skylab-origin pipelines

> **Important:** passing womtool validation is **not sufficient** to satisfy cascading bump requirements. Womtool only checks task call signatures — it does not enforce changelog updates. When CI reports "X.changelog.md has not been changed and needs to be updated", that is a cascading bump violation: patch-bump the WDL version, add a changelog entry (even "no functional impact"), and update `pipeline_versions.txt`.

## WDL Style

- **Version declaration:** WDL files declare `version 1.0`.
- **Formatting:** 2-space indentation; blank lines to separate logical sections; no strict line-length limit.
- **Naming:** tasks and call aliases use `UpperCamelCase`; variables use `lowercase_underscore` (Python style).
- **Workflow input block order:** required inputs first, optional inputs with defaults second, runtime-configuration parameters last.
- **Task section order:** input → command → output → runtime. In `command` blocks, put one input argument per line for clarity.
- **`meta { allowNestedInputs: true }`** — include for Terra compatibility.
- **GPU runtime keys:** use camelCase `gpuType`, `gpuCount`, and `nvidiaDriverVersion` in `runtime` blocks — universally portable across Cromwell/Terra/GCP. The snake_case aliases (`hardware_gpu_type`, `nvidia_driver_version`) are **not** portable.
- **WARP-tools vs inline WDL:** package complex or reusable processing scripts (Python, R, specialized deps, data transforms) into [warp-tools](https://github.com/broadinstitute/warp-tools) Docker images. Use inline WDL only for simple file operations, basic string manipulation, parameter validation, and straightforward conditional logic.
- **Inline Python in `command`:** use the unquoted `<<CODE` / `CODE` heredoc delimiter when invoking inline Python with `python3` — the established WARP convention. Do **not** use quoted delimiters (e.g. `<<'PYEOF'`). Use `~{wdl_var}` for WDL interpolation (resolved before bash runs); avoid f-strings combined with `~{}` interpolation since both use curly braces and conflict (prefer `'~{prefix}_output.txt'` over `f'~{prefix}_output.txt'`).

## Multi-cloud Support

Pipelines support GCP/Azure via a `cloud_provider` input. The pattern varies by pipeline origin:

**Broad-origin pipelines** (DNA-seq, arrays, genotyping) maintain docker registries and validate `cloud_provider` at runtime:

```wdl
String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
```

**Error handling:** use the `ErrorWithMessage` task for input validation. For broad-origin pipelines validating `cloud_provider`:

```wdl
if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
  call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
    input: message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
  }
}
```

## Pipeline Documentation Layout

Expect a pipeline's documentation to be split across two places. This isn't an enforced policy — it's just how pipelines evolve: a README starts in the pipeline directory, and the full write-up later lands on the docs site.

- **In-repo `README.md`** (`pipelines/wdl/<name>/README.md`) — usually slim and mostly a pointer to the docs-site page: a one-sentence summary, a link to the full docs page on the WARP site, a minimal "Running the pipeline" snippet, required inputs at a glance, and a link to `<Pipeline>.changelog.md`. Keep it that way — don't mirror the full docs page into it.
- **Docs-site page** (`website/docs/Pipelines/<Name>_Pipeline/README.md`) — the full documentation. Required Docusaurus frontmatter:

  ```yaml
  ---
  sidebar_position: 1
  slug: /Pipelines/<Name>_Pipeline/README
  ---
  ```

  Plus a sibling `_category_.json` when creating a new category.

When you change a pipeline's interface, update the docs-site page so its inputs/outputs stay accurate, and make sure the slim in-repo README still points correctly.

**Validate docs changes** with a full build — it catches broken links and bad frontmatter:

```bash
yarn --cwd=website install   # first time only
yarn --cwd=website build
```

`yarn --cwd=website start` is fine for previewing but does **not** fail on broken links — always use `build` before marking docs work complete.

**Cross-page links:** prefer relative paths to other pipeline pages (e.g. `[Multiome](../Multiome_Pipeline/README.md)`). When the target lacks a docs page, link to the GitHub source rather than a non-resolving relative path.

For Docusaurus markdown conventions follow [doc_style.md](website/docs/contribution/contribute_to_warp_docs/doc_style.md); for site build/serve commands follow [docsite_maintenance.md](website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md).

## CI / GitHub Workflows

The branch model (develop/staging/master) and how the smart test system selects which pipelines to test are documented in [TestingPipelines.md](website/docs/About_WARP/TestingPipelines.md). Operational specifics:

- **Path-based triggers** — GitHub Actions workflows trigger on the paths a pipeline depends on:

  ```yaml
  paths:
    - 'pipelines/wdl/optimus/**'
    - 'tasks/wdl/StarAlign.wdl'
    - 'tasks/wdl/Metrics.wdl'
    - 'tasks/wdl/Utilities.wdl'
  ```

  When adding a regression test input, wire its path into the relevant workflow so CI picks it up.
- **Test WDL naming & contract** — see [*Registering a CI test for a pipeline*](#registering-a-ci-test-for-a-pipeline-plumbing--scientific) above for the full `Test<Pipeline>.wdl` contract (naming, imports, `GetValidationInputs`, `update_truth`).
- **Release builds** — `./scripts/build_pipeline_release.sh -w pipeline.wdl -v version -o output_dir`.

## Git Merge Conflict Resolution

When resolving merge conflicts in WARP:

1. **WDL files first** — resolve WDL pipeline conflicts before configuration files.
2. **Version precedence** — use the higher version number from develop, then increment by 0.0.1.
3. **Import path updates** — check and fix import statements after resolving conflicts.
4. **Validation required** — run womtool validation on all modified WDL files.
5. **`pipeline_versions.txt` last** — resolve it last, as it's the master tracking file.

```bash
# 1. Resolve all .wdl files first
# 2. Handle .dockstore.yml
# 3. Process changelog files
# 4. Resolve pipeline_versions.txt last
# 5. Validate all WDL files with womtool
# 6. Commit and push
```

**After merging develop**, grep each touched WDL for duplicate `pipeline_version` declarations, which merge commits can silently introduce:

```bash
grep -n 'pipeline_version' pipelines/wdl/*/*.wdl | grep -v '\.changelog' | awk -F: '{print $1}' | sort | uniq -d
```

Keep the higher version and delete the duplicate line.

### Post-Merge Validation Checklist

- [ ] No legacy `tasks/broad/` or `tasks/skylab/` import paths remain
- [ ] No duplicate `pipeline_version` declarations in any WDL
- [ ] All WDL files pass womtool validation (run the loop above)
- [ ] Pipeline versions incremented and updated in `pipeline_versions.txt`
- [ ] Changelogs updated with current date and description
- [ ] `.dockstore.yml` paths reflect current structure

## Key Files and Directories

- `pipeline_versions.txt` — master version tracking
- `.dockstore.yml` — Dockstore registry configuration
- `.github/workflows/` — CI/CD workflows
- `wreleaser/` — CLI tool for querying releases
- `scripts/common.sh` — shared build functions
- `scripts/validate_wdls.sh` — convenience womtool wrapper
- `verification/test-wdls/scripts/` — test generation utilities

## Agent-Specific Notes

### `input_id` output prefix pattern

When a pipeline emits per-sample artifacts, accept a `String input_id` workflow input and prefix every output filename with `~{input_id}_`. This matches the convention used by scANVI, Optimus, Multiome, and other Skylab-origin pipelines.

### Optimus chemistry and reference assets

When adding support for a new 10x chemistry in Optimus:

- Use `tenx_chemistry_version` (Integer, e.g. `4`) for the major chemistry version and `tenx_chemistry_subversion` (optional String, e.g. `"v4_TRU"`) for whitelist variant selection within that version.
- Optimus whitelist files live at `gs://gcp-public-data--broad-references/optimus_whitelists/`. Do **not** use the old `RNA/resources/` path (which contained a `febrary` typo and is incorrect).
- Validate that inputs required by a given chemistry (e.g. `i1_fastq` for v4) are enforced via `ErrorWithMessage` when the chemistry version demands them.

### Sizing Optimus test data

Shrink Plumbing/Scientific inputs by subsetting reads to a genomic **region** (e.g. one chromosome), which keeps the called-cell count high. Do **not** subsample down to a handful of **barcodes**: the library-level doublet step sets its KNN neighbor count to `k = round(min(100, n_called_cells × 0.01))`, so below ~50 cells `k` becomes `0` and the run dies with sklearn's `n_neighbors ... Got 0` error. Keep ≥~500 called cells. See the [Minimum cell count warning](website/docs/Pipelines/Optimus_Pipeline/README.md).

### Jupyter notebooks

Never hand-write `.ipynb` JSON. Build notebooks programmatically with `nbformat` (markdown + code cells as plain strings), then **execute** them with `jupyter nbconvert --to notebook --execute --inplace …` so the delivered notebook has its figures and tables rendered inline. If the (often containerized) environment lacks `nbconvert`/`ipykernel`, install them there first — a scientific report container with `scanpy`/`matplotlib` usually has network for `pip`. **Deliver an *executed* notebook or stop and ask** — never ship an unexecuted notebook or one assembled by hand.

### When you discover a recurring agent mistake or lack of adherence

Record it here, under *Agent-Specific Notes* — this file is the only home for agent guidance. Keep entries short and link out to the canonical reference rather than duplicating it. **Before adding a new note, check whether the information is already covered above or in a linked doc** — prefer a one-line reference over a duplicate section.

If the user requests updating the AGENTS.md file with what you have learned, review this file and propose concise additions for future agents — see the reflection prompt in [AGENT_reflections.md](AGENT_reflections.md).
