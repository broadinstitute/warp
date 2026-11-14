# WDL authoring: agent guide

Read this for WDL work beyond a small edit. The general rules (validation, versioning, core style) are in [AGENTS.md](../AGENTS.md); the human-facing rationale is in [BestPractices.md](../website/docs/About_WARP/BestPractices.md).

## Validating a task and its importers

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

## Known cascade chains

A starting point only; `scripts/validate_release.sh` is authoritative (see [Cascading version bumps](../AGENTS.md#cascading-version-bumps)).

- `tasks/wdl/CheckInputs.wdl` → **Optimus** (`checkOptimusInput`) + its wrappers **Multiome**, **PairedTag**, **SlideTags**; **MultiSampleSmartSeq2SingleNucleus** (`checkInputArrays`); **SlideSeq** imports it but calls nothing (→ "no functional impact").
- `tasks/wdl/StarAlign.wdl` → **Optimus**, **SlideSeq**, **Multiome**, **PairedTag**, **SlideTags**, **MultiSampleSmartSeq2SingleNucleus**.
- `tasks/wdl/Metrics.wdl` → most Skylab-origin pipelines.

## Documented defaults

If a comment, `parameter_meta`, or README says an optional input "defaults to X when unspecified", make that default declarative (`select_first([input, "X"])` or a defaulted declaration) and verify the consuming code's unset path actually yields X. A `String?` interpolates to an empty string in bash, so a downstream `if/else` can silently encode the *opposite* default (as happened with Optimus `tenx_chemistry_subversion`); womtool checks types, not this, and tests that pin explicit values never exercise the default. Treat every prose "default" as a claim to verify against the code.

## Unused imports

womtool does **not** flag an `import "..." as NS` whose namespace is never used, so dead imports accumulate (they make a reader chase a dependency that isn't there). Scan for them:

```bash
scripts/find_unused_wdl_imports.sh
```

It skips `structs/` imports — WDL references struct types by bare name, so those legitimately have no `NS.` usage.

**Removing a dead import from a pipeline or shared `tasks/wdl/` WDL is still a WDL change**, so it forces a patch bump + changelog (+ cascade for shared tasks) — see [Cascading version bumps](../AGENTS.md#cascading-version-bumps). Dead imports in `verification/**` are free to drop. So don't do a repo-wide sweep in one PR: clean the `verification/**` ones anytime, and fold each pipeline/task import removal into that pipeline's **next** real change so it shares an already-required bump. Known debt as of this writing (run the scan for the current list): `atac.wdl` (`Merge`), `UltimaGenomicsWholeGenomeGermline.wdl` (`InternalTasks`, `QC`, `UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates`, `UltimaGenomicsWholeGenomeGermlineQC`), `UltimaGenomicsWholeGenomeCramOnly.wdl` (`VariantDiscoverTasks`), `SlideSeq.wdl` (`OptimusInputChecks`), `tasks/wdl/SplitLargeReadGroup.wdl` (`Utils`).

**Future CI wiring (not done — deliberately):** the test system is currently flaky, so this stays a manual scan rather than a gate. When CI is trusted, add it as an **advisory** (non-blocking, `continue-on-error`) job first so it reports without failing PRs; promote to a blocking check only once the known debt above is cleared and the signal is proven quiet.

## File naming

Never bake a lifecycle qualifier — `Updated`, `New`, `Old`, `Final`, `V2`, `Deprecated`, ... — into a WDL filename or its `workflow` name; that's changelog language, not identity. It shows up when a WDL is replaced but the replacement keeps a disambiguating suffix instead of being renamed once the original is deleted. Fix: delete the dead file and rename the replacement in the same change, updating every reference — the `import` path, the `as` alias, any call-site namespace, and CI `paths:` filters (see [Stale CI path filters](ci-tests.md#stale-ci-path-filters)). Example: `verification/VerifyCramToUnmappedBamsUpdated.wdl` outlived the `VerifyCramToUnmappedBams.wdl` it had replaced; a stale-migration artifact caught and fixed in #1913.

## GPU tasks

- **GPU runtime keys:** use camelCase `gpuType`, `gpuCount`, and `nvidiaDriverVersion` in `runtime` blocks — universally portable across Cromwell/Terra/GCP. The snake_case aliases (`hardware_gpu_type`, `nvidia_driver_version`) are **not** portable.
- **GPU-optional tasks need TWO task definitions.** Cromwell rejects `gpuCount: 0` ("Expecting gpuCount runtime attribute value greater than 0") *and* rejects optional-typed GPU attributes (`gpuCount: Int?` → "Expected positive Int but got Int? null"; `gpuType: String?` → "String value required but got String?"), and WDL 1.0 cannot conditionally omit a `runtime` key. So a single task is either always-GPU or never-GPU. To offer both GPU and CPU-only runs, write two tasks — identical except one has the GPU `runtime` attributes and the other omits them — and route with a workflow `if (gpu_count > 0)` / `if (gpu_count == 0)`, gathering outputs with `select_first`. Keep the container **device-agnostic** (auto-detect the accelerator) so both tasks run the same command and only VM *provisioning* differs. Example: `pipelines/wdl/scanvi/scANVI.wdl` (`MultiomeLabelTransfer` + `MultiomeLabelTransferCpu`).

## Preemptibility

Tasks wire a `preemptible_tries` input to `preemptible:` in `runtime`. Keep it preemptible for short/retryable steps, but pass `preemptible_tries = 0` (non-preemptible) for long single-threaded steps whose runtime can exceed the preemption window — e.g. a large-BAM gather/merge — because a preemption wastes the whole run. Signature in a failed test: the task hits repeated `RetryableFailure` and can end on a VM-teardown `exit 125` ("job was stopped before the command finished") while the tool itself had actually completed. Precedent: `Qc.wdl`, `JointGenotypingTasks.wdl`, `UnmappedBamToAlignedBam.wdl` (`GatherBamFiles`).

## `command` block examples

For bash scripts, use the `<<<` / `>>>` style:

```wdl
command <<<
  set -e
  echo "Hello, ~{name}!"
>>>
```

For inline Python, use the `<<CODE` / `CODE` style. This is the established convention for tasks that primarily wrap a Python script:

```wdl
command <<CODE
  set -e
  python3 -c '
  # Python code here.
  # Use ~{wdl_var} for WDL interpolation.
  '
CODE
```

## Multi-cloud support

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

## Optimus chemistry and reference assets

When adding support for a new 10x chemistry in Optimus:

- Use `tenx_chemistry_version` (Integer, e.g. `4`) for the major chemistry version and `tenx_chemistry_subversion` (optional String, e.g. `"v4_TRU"`) for whitelist variant selection within that version.
- Optimus whitelist files live at `gs://gcp-public-data--broad-references/optimus_whitelists/`. Do **not** use the old `RNA/resources/` path (which contained a `febrary` typo and is incorrect).
- Validate that inputs required by a given chemistry (e.g. `i1_fastq` for v4) are enforced via `ErrorWithMessage` when the chemistry version demands them.
