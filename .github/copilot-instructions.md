# WARP AI Coding Guidelines

## Architecture Overview

WARP is a collection of cloud-optimized WDL (Workflow Description Language) pipelines for biological data processing. The repository follows a structured pattern with a **flattened directory structure**:

- **`pipelines/wdl/`**: All workflow definitions in a unified directory structure
- **`tasks/wdl/`**: All reusable WDL tasks in a unified directory structure
- **`structs/`**: WDL struct definitions for type safety
- **`verification/`**: Test workflows that validate pipeline outputs
- **`scripts/`**: Build and validation automation

```
pipelines/wdl/
├── arrays/
├── atac/
├── dna_seq/
├── genotyping/
├── multiome/
├── optimus/
├── peak_calling/
├── reprocessing/
├── rna_seq/
└── ... (other pipeline categories)

tasks/wdl/
├── Alignment.wdl
├── BamProcessing.wdl
├── FastqProcessing.wdl
├── StarAlign.wdl
├── Utilities.wdl
└── ... (all task files)
```

> **History:** Content was previously split between `pipelines/broad/` (DNA-seq, arrays, reprocessing) and `pipelines/skylab/` (single-cell RNA-seq, ATAC-seq, multiome). Both are now unified under `pipelines/wdl/` and `tasks/wdl/`. Docker images are maintained separately in [warp-tools](https://github.com/broadinstitute/warp-tools).

## WDL Development Patterns

### Import Structure
All WDL files use relative imports:
```wdl
import "../../../tasks/wdl/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../structs/dna_seq/DNASeqStructs.wdl"
```

Verify no legacy paths remain: the following command should return nothing.
```bash
grep -rn 'tasks/broad\|tasks/skylab\|pipelines/broad\|pipelines/skylab' pipelines tasks verification
```

### Sub-workflow Input Contract
A WDL file defines a **single input contract for all callers**. You cannot expose an input to one workflow but hide it from another that imports the same WDL. To remove a parameter from one consumer, you must remove it from the shared task/sub-workflow **and** from every caller.

### Pipeline Versioning
Every pipeline has a semantic version string. Bump the patch version **once per branch/PR**:
```wdl
workflow MyPipeline {
  String pipeline_version = "1.2.4"  # Increment for changes
}
```

**If the top changelog entry on the current branch is already an unreleased bump**, append bullets to that entry rather than creating a new version. Only `pipeline_versions.txt` and the WDL `pipeline_version` string need to stay in sync with the top of the changelog.

### Changelog Requirements
Each pipeline has a `<PipelineName>.changelog.md` file. For each new version:
1. Add entry at the top with format: `# 1.2.5\nYYYY-MM-DD (Date of Last Commit)\n\n* Description of changes`
2. Use past tense, active voice: "Updated docker image" not "Docker image was updated"
3. **Date Format**: Always use YYYY-MM-DD matching the commit date. The literal suffix `(Date of Last Commit)` is part of the format.

**Example Changelog Entry:**
```markdown
# 1.2.5
2026-04-28 (Date of Last Commit)

* Updated docker image to latest version
* Fixed validation issue with input parameters

# 1.2.4
2026-03-15 (Date of Last Commit)

* Added new optional input
```

### Cascading Version Bumps
When a shared task changes, every workflow that imports it may need a version bump. Common dependency chains:

- `tasks/wdl/StarAlign.wdl` → **Optimus**, **SlideSeq**, **Multiome** (via Optimus), **PairedTag** (via Optimus), **SlideTags** (via Optimus), **MultiSampleSmartSeq2SingleNucleus**
- `tasks/wdl/CheckInputs.wdl` → **Optimus** (via `checkOptimusInput`) and its wrappers (**Multiome**, **PairedTag**, **SlideTags**); also **MultiSampleSmartSeq2SingleNucleus** (via `checkInputArrays`, a separate task in the same file). Note: SlideSeq imports this file but never calls any of its tasks — changes to CheckInputs.wdl task signatures do not functionally affect SlideSeq through this import.
- `tasks/wdl/Metrics.wdl` → Most Skylab-origin pipelines

For each dependent pipeline, either bump the version with a changelog note **or** add a "no functional impact" note explaining why the change is benign for that pipeline.

> **Important:** passing womtool validation is **not sufficient** to satisfy cascading bump requirements. Womtool only checks task call signatures — it does not enforce changelog updates. When CI reports "X.changelog.md has not been changed and needs to be updated", that is a cascading bump violation: patch-bump the WDL version, add a changelog entry (even "no functional impact"), and update `pipeline_versions.txt`.

### Test Inputs
Test JSON inputs live at two locations per pipeline:
- `pipelines/wdl/<name>/test_inputs/{Plumbing,Scientific}/` — CI test inputs (Plumbing = fast/cheap, Scientific = end-to-end accuracy)
- `pipelines/wdl/<name>/example_inputs/` — example inputs for documentation

When adding a new test input for a regression case, also wire it into the relevant GitHub Actions workflow under `.github/workflows/` so CI picks it up.

### Testing Pattern
Test workflows follow strict naming: `Test<PipelineName>.wdl` in `verification/test-wdls/`
- Import main pipeline and verification workflow
- Use `GetValidationInputs` tasks to compare test vs truth outputs
- Include `update_truth` boolean for refreshing golden files

### Release Process
Use build scripts for releases:
```bash
./scripts/build_pipeline_release.sh -w pipeline.wdl -v version -o output_dir
```

## Critical Workflows

### WDL Validation (MANDATORY)
The womtool JAR is expected at a user-specific location such as `~/womtool-90.jar`. The appropriate version of Java must also be available. A convenience wrapper also exists at `scripts/validate_wdls.sh`.

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

**Critical Notes:**
- Womtool validation is REQUIRED for all WDL changes
- Validation must pass for all modified pipeline files **and their dependents**
- Import path errors are the most common validation failures

### Pipeline Version Management
1. **Version Bump Rule**: Increment patch version (e.g., 1.2.3 → 1.2.4) at *most* once per branch/PR per merge from develop to staging.
2. **Pipeline-Specific Updates**: Update version string in the WDL file
3. **Master Tracking**: Update `pipeline_versions.txt` with new version and current date
4. **Linked Versions**: Some pipelines reference sub-workflow versions (e.g., ImputationBeagle references ArrayImputationQuotaConsumed version)

## Structural Conventions

### Multi-cloud Support
Pipelines support GCP/Azure via a `cloud_provider` input. The pattern varies by pipeline origin:

**Broad-origin pipelines** (DNA-seq, arrays, genotyping) maintain dual docker registries and validate `cloud_provider` at runtime:
```wdl
String gatk_docker_gcp = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
String gatk_docker_azure = "terrapublic.azurecr.io/gatk:4.6.1.0"
String gatk_docker = if cloud_provider == "gcp" then gatk_docker_gcp else gatk_docker_azure
```

**Skylab-origin pipelines** (Optimus, Multiome, ATAC-seq, PairedTag, SlideTags, SlideSeq, RNAWithUMIs, MultiSampleSmartSeq2) default `cloud_provider = "gcp"` and do **not** maintain Azure docker registries. Do not add Azure docker alternatives to these pipelines.

### Error Handling
Use `ErrorWithMessage` task for input validation. For broad-origin pipelines validating `cloud_provider`:
```wdl
if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
  call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
    input: message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
  }
}
```

### Meta Parameters
Include for Terra compatibility:
```wdl
meta {
  allowNestedInputs: true
}
```

### WDL Style
Follow [WARP_WDL_Style_Guide.md](../WARP_WDL_Style_Guide.md) for formatting, naming, `parameter_meta`, and the `<<CODE` heredoc convention for inline Python in `command` blocks (quoted terminators like `<<'PYEOF'` suppress `~{}` interpolation and are not used here).

### GPU Runtime Keys
Use camelCase `gpuType`, `gpuCount`, and `nvidiaDriverVersion` in `runtime` blocks — universally portable across Cromwell/Terra/GCP. The snake_case aliases (`hardware_gpu_type`, `nvidia_driver_version`) are not portable.

## Documentation

Two-tier layout for every pipeline:

- **In-repo `README.md`** (`pipelines/wdl/<name>/README.md`) — slim: one-sentence summary, link to the full docs page on the WARP site, minimal "Running the pipeline" snippet, required inputs at a glance, link to the pipeline changelog. Do not mirror the full docs page here.
- **Docs-site page** (`website/docs/Pipelines/<Name>_Pipeline/README.md`) — full documentation. Required Docusaurus frontmatter:
  ```yaml
  ---
  sidebar_position: 1
  slug: /Pipelines/<Name>_Pipeline/README
  ---
  ```

For Docusaurus markdown conventions (admonitions, tabs, code-block highlighting, tables, cross-refs) follow [website/docs/contribution/contribute_to_warp_docs/doc_style.md](../website/docs/contribution/contribute_to_warp_docs/doc_style.md). For site build/serve commands follow [website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md](../website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md).

**Validate docs changes** with a full build (catches broken links and bad frontmatter):
```bash
yarn --cwd=website install   # first time only
yarn --cwd=website build
```
`yarn --cwd=website start` is fine for previewing but does not fail on broken links.

**Cross-page links**: prefer relative paths to other pipeline pages (e.g. `[Multiome](../Multiome_Pipeline/README.md)`). When the target lacks a docs page, link to the GitHub source rather than a non-resolving relative path.

**Stale example/test inputs**: when renaming a workflow or removing/renaming inputs, audit `pipelines/wdl/<name>/example_inputs/*.json` and `test_inputs/**/*.json` — these break silently because they are not validated by womtool.

## GitHub Workflows and CI/CD

### Workflow Path Triggers
GitHub Actions workflows use path-based triggers:
```yaml
paths:
  - 'pipelines/wdl/optimus/**'
  - 'tasks/wdl/StarAlign.wdl'
  - 'tasks/wdl/Metrics.wdl'
  - 'tasks/wdl/Utilities.wdl'
```

### Testing System
- Tests are triggered when files change in relevant `pipelines/wdl/` or `tasks/wdl/` paths
- Each pipeline has a dedicated test workflow file
- Verification WDLs in `verification/` directory validate outputs

## Key Files and Directories

- `pipeline_versions.txt`: Master version tracking
- `changelog_style.md`: Detailed changelog formatting rules
- `.dockstore.yml`: Dockstore registry configuration
- `.github/workflows/`: CI/CD workflows
- `wreleaser/`: CLI tool for querying releases
- `scripts/common.sh`: Shared build functions
- `scripts/validate_wdls.sh`: Convenience womtool wrapper
- `verification/test-wdls/scripts/`: Test generation utilities

## Git Merge Conflict Resolution

When resolving merge conflicts in WARP:

1. **WDL Files First**: Always resolve WDL pipeline conflicts before configuration files
2. **Version Precedence**: Use higher version numbers from develop branch, then increment by 0.0.1
3. **Import Path Updates**: Check and fix import statements after resolving conflicts
4. **Validation Required**: Run womtool validation on all modified WDL files
5. **Pipeline Versions Last**: Resolve `pipeline_versions.txt` conflicts last as it's the master tracking file

**Conflict Resolution Order:**
```bash
# 1. Resolve all .wdl files first
# 2. Handle .dockstore.yml
# 3. Process changelog files
# 4. Resolve pipeline_versions.txt last
# 5. Validate all WDL files with womtool
# 6. Commit and push
```

**After merging develop**, grep each touched WDL for duplicate `pipeline_version` declarations, which can be silently introduced by merge commits:
```bash
grep -n 'pipeline_version' pipelines/wdl/*/*.wdl | grep -v '\.changelog' | awk -F: '{print $1}' | sort | uniq -d
```
Keep the higher version and delete the duplicate line.

### Post-Merge Validation Checklist
After completing a merge:
- [ ] No legacy `tasks/broad/` or `tasks/skylab/` import paths remain
- [ ] No duplicate `pipeline_version` declarations in any WDL
- [ ] All WDL files pass womtool validation (run the loop above)
- [ ] Pipeline versions incremented and updated in `pipeline_versions.txt`
- [ ] Changelogs updated with current date and description
- [ ] `.dockstore.yml` paths reflect current structure
