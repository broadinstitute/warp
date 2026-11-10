# Agent Instructions

Single entry point for all agent-driven and AI-coding work in this repository. Start here. This file holds the rules most tasks need; task-specific guidance lives in the guides under [agents/](agents/), listed below. For the *rules* and *rationale* behind both, follow the links to the canonical human-facing docs — those are the source of truth, so do not restate them here.

> **Opening, updating, or merging pull requests is never an agent's job.** Do the work on a branch and push it when asked; the human opens and manages the PR. Also confirm a branch still exists before pushing to it — pushing to a deleted (e.g. post-merge) branch silently recreates it.

## Task guides

Read the matching guide before you start:

| Task | Guide |
| --- | --- |
| Debugging a failed Terra/WDL pipeline test (submission URL → the exact failing task/metric, via the `fiss-mcp` tools) | [agents/fiss-mcp-runbook.md](agents/fiss-mcp-runbook.md) |
| Adding or changing a Plumbing/Scientific test, test inputs, a `Verify*` compare task, or a `.github/workflows/test_*.yml` | [agents/ci-tests.md](agents/ci-tests.md) |
| WDL work beyond a small edit: GPU tasks, preemptibility, multi-cloud, renaming files, dead imports, optional-input defaults, Optimus chemistry | [agents/wdl-authoring.md](agents/wdl-authoring.md) |
| Pipeline docs pages and the Docusaurus site | [agents/pipeline-docs.md](agents/pipeline-docs.md) |

## Authoritative References

| Topic | Document |
| --- | --- |
| Changelog format and language (entry structure, active voice, sample entries) | [changelog_style.md](website/docs/contribution/contribute_to_warp/changelog_style.md) |
| Semantic-versioning classification (major/minor/patch) and release packaging | [VersionAndReleasePipelines.md](website/docs/About_WARP/VersionAndReleasePipelines.md) |
| Pipeline best practices (automation, testability, portability, scaling, maintainability) | [BestPractices.md](website/docs/About_WARP/BestPractices.md) |
| Requirements a released pipeline must meet | [PipelineRequirements.md](website/docs/About_WARP/PipelineRequirements.md) |
| Testing: branch model (develop/staging/master) and how the smart test system selects pipelines | [TestingPipelines.md](website/docs/About_WARP/TestingPipelines.md) |
| Docusaurus markdown style (admonitions, tables, tabs, code-block highlighting, cross-refs) | [doc_style.md](website/docs/contribution/contribute_to_warp_docs/doc_style.md) |
| How to build / serve / validate the docs site | [docsite_maintenance.md](website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md) |

## Quick Checklist Before Completing Any WARP Task

1. **Validate all modified WDLs** — and every WDL that imports them (transitively). See [WDL Validation (MANDATORY)](#wdl-validation-mandatory).
2. **Changelog and versioning** — bump versions and update changelogs per [Changelog and Versioning](#changelog-and-versioning), including the cascade to every pipeline that imports a changed shared task. Validation passing is not sufficient.
3. **After merging develop**, check for duplicate `pipeline_version` lines — keep the higher one. See [Git](#git).
4. **Sub-workflow contract** — a WDL file defines a single input contract for all callers: you cannot expose an input to one workflow but hide it from another that imports the same WDL. Removing an input from a shared task/sub-workflow requires removing it from every caller.
5. **Stale example/test inputs** — when you rename a workflow or add, remove or rename inputs, audit `example_inputs/*.json` and `test_inputs/**/*.json` and thread the change through the `Test<Name>.wdl` wrapper; they break silently because womtool does not check them. See [ci-tests.md → Test inputs](agents/ci-tests.md#test-inputs).
6. **Touching a pipeline's interface** — also update its docs page under `website/docs/Pipelines/<Name>_Pipeline/README.md`, link version information from the changelog instead of hardcoding it, and run `yarn --cwd=website build` to catch broken links. See [pipeline-docs.md](agents/pipeline-docs.md).
7. **Documented defaults must match the code** — if a comment, `parameter_meta`, or README says an optional input "defaults to X when unspecified", make that default declarative (`select_first([input, "X"])` or a defaulted declaration) and verify the unset path actually yields X. See [wdl-authoring.md → Documented defaults](agents/wdl-authoring.md#documented-defaults).

## Repository Structure

WARP is a collection of cloud-optimized WDL pipelines for biological data processing, in a **flattened directory structure**:

- **`pipelines/wdl/`** — all workflow definitions, by category (`arrays/`, `dna_seq/`, `multiome/`, `optimus/`, …). A pipeline's directory holds its WDL, `<Pipeline>.changelog.md`, `test_inputs/` and `example_inputs/`. `pipelines/nextflow/` is a placeholder.
- **`tasks/wdl/`** — all reusable WDL tasks; **`structs/`** — WDL struct definitions.
- **`verification/`** — `Verify<Pipeline>.wdl` compare workflows; `verification/test-wdls/` holds the `Test<Pipeline>.wdl` CI wrappers.
- **`scripts/`** — build and validation automation (`scripts/common.sh` holds shared build functions); **`wreleaser/`** — CLI for querying releases.
- **`website/`** — Docusaurus documentation site; **`all_of_us/`** — All of Us Research Program workflows.
- **`.dockstore.yml`** — Dockstore registry; **`.github/workflows/`** — CI; **`pipeline_versions.txt`** — generated version tracking (never edit by hand, see below).

> **History:** content used to be split between `pipelines/broad/` and `pipelines/skylab/` (and `tasks/broad/`, `tasks/skylab/`); any remnants are ones you should offer to clean up. Docker images are maintained separately in [warp-tools](https://github.com/broadinstitute/warp-tools).

## WDL Validation (MANDATORY)

Womtool validation is REQUIRED for all WDL changes, and must pass for every modified file **and its dependents** — a task can validate in isolation but break a caller's argument list. The womtool JAR is expected at a user-specific location such as `~/womtool-90.jar` (with a compatible Java available; if missing: `curl -sSLo ~/womtool-90.jar https://github.com/broadinstitute/cromwell/releases/download/90/womtool-90.jar`); a convenience wrapper exists at `scripts/validate_wdls.sh`.

```bash
java -jar ~/womtool-90.jar validate path/to/pipeline.wdl
```

For a loop that validates a task and all its importers, see [wdl-authoring.md](agents/wdl-authoring.md#validating-a-task-and-its-importers). womtool checks types and call signatures only — not test JSONs, changelogs, unused imports, CI YAML, or `.dockstore.yml`.

## Changelog and Versioning

Every pipeline carries a `String pipeline_version = "major.minor.patch"` and a cumulative `<PipelineName>.changelog.md`. For entry **format and language** follow the [changelog style guide](website/docs/contribution/contribute_to_warp/changelog_style.md); for what makes a change **major / minor / patch** follow [VersionAndReleasePipelines.md](website/docs/About_WARP/VersionAndReleasePipelines.md). The agent-operational rules below are not in those docs:

1. **Version bumps are once per branch/PR.** If the top changelog entry on the current branch already has an unreleased version bump (i.e. its version is higher than what is in `develop`), append new bullet points to that entry — do not create a new version section.
2. **`pipeline_versions.txt`** is automatically generated by a GitHub Action (`.github/workflows/update_versions_file.yml`) on pull requests to `develop`. **Do not edit this file manually.** The action **commits to your PR branch** after any push that changes a version or changelog date, so `git pull` before the next commit/push or the push is rejected as non-fast-forward. The action derives its content from the version and date at the top of each pipeline's `.changelog.md`.
3. **Linked versions** — some pipelines reference a sub-workflow's version (e.g. ImputationBeagle references the ArrayImputationQuotaConsumed version); update both together.
4. **What requires a version bump.** A change to a pipeline WDL (`pipelines/wdl/**`) or a shared `tasks/wdl/` task it imports requires at least a patch `+0.0.1` bump on every affected pipeline, classified per pipeline from the guide (e.g. adding an [optional input](website/docs/About_WARP/VersionAndReleasePipelines.md#patch-micro-changes), one with a default or a `?` type, is **patch** when leaving it out keeps the outputs unchanged; a change that may produce slightly different outputs, such as a rebuilt Docker image, is **minor**), **and** a changelog entry. A change confined to `verification/**` — the `Test<Pipeline>.wdl` wrappers and the `Verify*` / `CompareMetrics`-style compare tasks — requires **no** `pipeline_version` bump and **no** changelog, because it ships no pipeline artifact. (Still validate it with womtool.)

### Cascading version bumps

Editing a shared task (`tasks/wdl/*.wdl`) forces a version bump on **every** pipeline that imports it, directly or transitively. For each such pipeline update **two** things: the WDL's `String pipeline_version`, and its `.changelog.md` (new entry — a plain "No functional impact" bullet is fine when the change can't affect it).

**Get the exact list from the tool — don't guess:** run `scripts/validate_release.sh -g origin/staging` (add `-i true` to also check changelogs). It prints `X.wdl has not been changed and needs updating` for every importer still missing a bump; fix and re-run until it says all are valid. This is exactly what the `WARP Validate Version` / `WARP Validate Changelog` CI checks run, so a local pass means those checks pass. Known import chains are listed in [wdl-authoring.md](agents/wdl-authoring.md#known-cascade-chains), but the script is authoritative.

> **Commit first.** The script diffs committed `HEAD` against the ref (`git diff HEAD <ref>`), so uncommitted working-tree edits are invisible — with changes unstaged it will falsely report all-valid. Commit, then run it. (It *does* trace task→pipeline imports transitively, so once committed a shared-task change correctly flags every importer.)

## WDL Style

- **Imports** are relative (`import "../../../tasks/wdl/UnmappedBamToAlignedBam.wdl" as ToBam`); import path errors are the most common validation failure.
- **Version and formatting:** `version 1.0`; 2-space indentation; blank lines between logical sections; no strict line-length limit.
- **Naming:** tasks and call aliases use `UpperCamelCase`; variables use `lowercase_underscore`. Never put a lifecycle qualifier (`Updated`, `New`, `V2`, …) in a WDL filename or `workflow` name — see [wdl-authoring.md](agents/wdl-authoring.md#file-naming).
- **Workflow input block order:** required inputs first, [optional inputs](website/docs/About_WARP/VersionAndReleasePipelines.md#patch-micro-changes) (a default or a `?` type) second, runtime-configuration parameters last.
- **Task section order:** input → command → output → runtime. In `command` blocks, put one input argument per line.
- **`meta { allowNestedInputs: true }`** — include for Terra compatibility.
- **`command` blocks:** use `<<<` / `>>>` for bash and `<<CODE` / `CODE` for tasks that wrap inline Python. Never quote the delimiter (`command <<'EOF'` disables WDL interpolation), and avoid Python f-strings with `~{}` interpolation, since both use curly braces (prefer `'~{prefix}_output.txt'`).
- **Input validation errors** use the `ErrorWithMessage` task (`tasks/wdl/Utilities.wdl`).
- **`input_id` output prefix:** a pipeline that emits per-sample artifacts accepts a `String input_id` and prefixes every output filename with `~{input_id}_` (as scANVI, Optimus and Multiome do).
- **WARP-tools vs inline WDL:** package complex or reusable processing scripts (Python, R, specialized deps, data transforms) into [warp-tools](https://github.com/broadinstitute/warp-tools) Docker images. Use inline WDL only for simple file operations, basic string manipulation, parameter validation, and straightforward conditional logic.

## Git

- **New branches track `origin/develop` by default.** `git checkout -b <branch> origin/develop` sets the upstream to `origin/develop`, so a plain `git push` refuses and `git pull` reports "Already up to date" against develop. Use `--no-track` (or push once with `git push -u origin <branch>`).
- **Merge conflicts:** resolve `.wdl` files first, then `.dockstore.yml`, then changelogs. Use the higher version number from develop, then increment by 0.0.1; fix import paths; update changelogs with the current date and description; validate every modified WDL with womtool.
- **After merging develop**, grep for duplicate `pipeline_version` declarations, which merge commits can silently introduce; keep the higher version and delete the duplicate line. Also check that no `tasks/broad/` or `tasks/skylab/` import paths came back and that `.dockstore.yml` paths match the current structure.

  ```bash
  grep -rc --include='*.wdl' 'String pipeline_version *=' pipelines/wdl | grep -v ':[01]$'
  ```

## Jupyter Notebooks

Never hand-write `.ipynb` JSON. Build notebooks programmatically with `nbformat` (markdown + code cells as plain strings), then **execute** them with `jupyter nbconvert --to notebook --execute --inplace …` so the delivered notebook has its figures and tables rendered inline. If the (often containerized) environment lacks `nbconvert`/`ipykernel`, install them there first — a scientific report container with `scanpy`/`matplotlib` usually has network for `pip`. **Deliver an *executed* notebook or stop and ask** — never ship an unexecuted notebook or one assembled by hand.

## Recording New Agent Guidance

When you discover a recurring agent mistake or lack of adherence, record it in the guide under [agents/](agents/) it belongs to, or in this file only if most tasks need it. Keep entries short and link out to the canonical reference rather than duplicating it. **Before adding a note, check whether it is already covered here, in a guide, or in a linked doc.**

If the user requests updating these files with what you have learned, follow the reflection prompt in [agents/reflections.md](agents/reflections.md).
