# WARP AI Coding Guidelines

These guidelines are for AI agents to interact with the WARP repo. They should also be a useful point of reference for humans too! Every time an agent interacts with WARP, please ask it to update these instructions to reflect on what it has learned, and to clarify and condense the instructions to reduce mistakes.

## Architecture Overview

WARP is a collection of cloud-optimized WDL (Workflow Description Language) pipelines for biological data processing. The repository follows a structured pattern with a **flattened directory structure**:

- **`pipelines/wdl/`**: All workflow definitions in a unified directory structure
- **`tasks/wdl/`**: All *reusable* WDL tasks in a unified directory structure  
- **`structs/`**: WDL struct definitions for type safety
- **`verification/`**: Test workflows that validate pipeline outputs
- **`scripts/`**: Build and validation automation

## Team History and Pipeline Organization

- **Historical Teams**: Content originally organized by Broad and Skylab teams
  - Broad team: DNA-seq, arrays, reprocessing 
  - Skylab team: Single-cell RNA-seq, ATAC-seq, multiome
- **Current Organization**: Unified under `pipelines/wdl/` and `tasks/wdl/` directories
- Docker images maintained separately in [warp-tools](https://github.com/broadinstitute/warp-tools)

### Flattened Directory Structure

As of the latest reorganization, the repository uses a flattened structure:

- **Old structure**: `pipelines/broad/` and `pipelines/skylab/` → **New**: `pipelines/wdl/`
- **Old structure**: `tasks/broad/` and `tasks/skylab/` → **New**: `tasks/wdl/`

This maintains the same directory depth while consolidating historical team-specific directories.

### Current Directory Structure

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
## Migration and Reorganization Notes

When working with import statements or file paths:

- All `pipelines/broad/` references should use `pipelines/wdl/`
- All `pipelines/skylab/` references should use `pipelines/wdl/`
- All `tasks/broad/` references should use `tasks/wdl/`
- All `tasks/skylab/` references should use `tasks/wdl/`
- Relative imports maintain the same depth but point to unified directories
- GitHub workflow path triggers have been updated accordingly

## WDL Development Patterns

### Import Structure

All WDL files use relative imports with this hierarchy:

```wdl
import "../../../tasks/wdl/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../structs/dna_seq/DNASeqStructs.wdl"
```

### Import Path Management

When modifying WDL files, ensure all import statements use the current structure:

```wdl
# Correct (current structure)
import "../../../tasks/wdl/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../tasks/wdl/Utilities.wdl" as utils

# Incorrect (legacy structure - will cause validation failures)
import "../../../tasks/broad/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../tasks/skylab/Utilities.wdl" as utils
```

**Automated Import Fix Script:**
For bulk import path updates across multiple files:

```bash
find pipelines/wdl -name "*.wdl" -exec sed -i 's|tasks/broad/|tasks/wdl/|g' {} \;
find pipelines/wdl -name "*.wdl" -exec sed -i 's|tasks/skylab/|tasks/wdl/|g' {} \;
```

Note that `-exec sed` commands are powerful and can modify many files at once. Use with caution and verify changes.

### Pipeline Versioning

Every pipeline has a semantic version string that MUST be incremented for any changes:

```wdl
workflow MyPipeline {
  String pipeline_version = "1.2.4"  # Increment for changes
}
```

Every pipeline modification requires version increment:

1. **Version Bump Rule**: Increment patch version (e.g., 1.2.3 → 1.2.4) for any change. Read the WARP semantic versioning policy to understand when to increment minor or major versions: `website/docs/About_WARP/VersionAndReleasePipelines.md`
2. **Pipeline-Specific Updates**: Update version string in the WDL file
3. **Master Tracking**: Update `pipeline_versions.txt` with new version and current date
4. **Linked Versions**: Some pipelines reference sub-workflow versions (e.g., ImputationBeagle references ArrayImputationQuotaConsumed version)

Example version update pattern:

```wdl
workflow MyPipeline {
  String pipeline_version = "1.2.5"  # Incremented from 1.2.4
}
```

### Changelog Requirements

Each pipeline has a `<PipelineName>.changelog.md` file. For every version bump:

1. Add entry at top with format: `# 1.2.5\n2025-MM-DD (Date of Last Commit)\n\n* Description of changes`
2. Use past tense, active voice: "Updated docker image" not "Docker image was updated"
3. **Standard Messages**: Use "Reorganize WDL pipelines into unified directory" for structural changes
4. **Date Format**: Always use YYYY-MM-DD format matching commit date

**Example Changelog Entry:**

```markdown
# 1.2.5
2025-09-19 (Date of Last Commit)

* Reorganize WDL pipelines into unified directory
* Updated docker image to latest version

# 1.2.4
2025-09-15 (Date of Last Commit)

* Fixed validation issue with input parameters
```

## Critical Workflows

### WDL Validation (MANDATORY - NO EXCEPTIONS)

**⚠️ CRITICAL REQUIREMENT: EVERY WDL file modification MUST be validated with womtool before committing ⚠️**

Always validate WDL syntax immediately after making ANY change to a WDL file. You need to know where womtool is located on your system. If you can't find it, ask your operator before proceeding. For example:

```bash
java -jar ~/womtool-90.jar validate path/to/pipeline.wdl
```

**Non-Negotiable Rules:**

- **EVERY WDL change** requires womtool validation - no exceptions, even for simple version updates. **NEVER commit** a WDL file without successful womtool validation
- **ALL modified WDL files** must pass validation, including dependencies
- **IMMEDIATE validation** - run womtool right after each WDL edit, not at the end
- Import path errors are the most common validation failures after directory reorganization
- If validation fails, check import paths first - they must use the current flattened structure

**Validation Workflow:**
1. Edit WDL file
2. **IMMEDIATELY** run `java -jar ~/womtool-90.jar validate <file.wdl>`
3. Fix any validation errors before proceeding
4. Only commit after successful validation

### Testing Pattern

Test workflows follow strict naming: `Test<PipelineName>.wdl` in `verification/test-wdls/`

- Import main pipeline and verification workflow
- Use `GetValidationInputs` tasks to compare test vs truth outputs
- Only include `update_truth` boolean for refreshing golden files if indicated to do so by a human

### Release Process

Use build scripts for releases:

```bash
./scripts/build_pipeline_release.sh -w pipeline.wdl -v version -o output_dir
```

## Structural Conventions

### Multi-cloud Support

Pipelines support GCP/Azure via `cloud_provider` input and conditional docker selection:

```wdl
String gatk_docker_gcp = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
String gatk_docker_azure = "terrapublic.azurecr.io/gatk:4.6.1.0"  
String gatk_docker = if cloud_provider == "gcp" then gatk_docker_gcp else gatk_docker_azure
```

### Error Handling

Use `ErrorWithMessage` task for input validation:

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

## GitHub Workflows and CI/CD

### Workflow Path Triggers

GitHub Actions workflows use path-based triggers for the flattened structure:

```yaml
paths:
  - 'pipelines/wdl/optimus/**'
  - 'tasks/wdl/StarAlign.wdl'
  - 'tasks/wdl/Metrics.wdl'
  - 'tasks/wdl/Utilities.wdl'
```

### Testing System

- Tests are triggered when files change in relevant `pipelines/wdl/` or pipeline associated `tasks/wdl/` paths
- Each pipeline has a dedicated test workflow file
- Verification WDLs in `verification/` directory validate outputs

## Key Files and Directories

### Configuration Files

- `pipeline_versions.txt`: Master version tracking
- `changelog_style.md`: Detailed changelog formatting rules  
- `.dockstore.yml`: Dockstore registry configuration (updated for flattened structure)
- `.github/workflows/`: CI/CD workflows (updated for new paths)

### Tools and Scripts

- `wreleaser/`: CLI tool for querying releases
- `scripts/common.sh`: Shared build functions
- `verification/test-wdls/scripts/`: Test generation utilities

### Git Merge Conflict Resolution

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

### Post-Merge Validation Checklist

After completing a merge:

- [ ] All import paths use `tasks/wdl/` structure
- [ ] All WDL files pass womtool validation
- [ ] Pipeline versions incremented and updated in `pipeline_versions.txt`
- [ ] **Version Consistency**: Verify WDL pipeline_version matches changelog version
- [ ] Changelogs updated with current date and description
- [ ] `.dockstore.yml` paths reflect current structure

### Version Consistency Validation

**Critical Issue**: WDL files and their corresponding changelog files can become out of sync during merges or manual edits. This creates confusion and violates WARP's version management principles.

**Validation Process**:
1. Always verify that `pipeline_version` in WDL matches the top entry in `<PipelineName>.changelog.md`
2. Use the higher version number when mismatches are found
3. If WDL version is higher: Add new changelog entry for that version
4. If changelog version is higher: Update WDL to match the changelog version
5. Never leave version mismatches unresolved

**Detection Command**:
```bash
# Compare WDL versions with changelog versions
for wdl in $(find pipelines/wdl -name "*.wdl"); do
  wdl_version=$(grep -Po 'pipeline_version = "\K[^"]*' "$wdl" 2>/dev/null)
  changelog="${wdl%.*}.changelog.md"
  if [[ -f "$changelog" ]]; then
    changelog_version=$(head -1 "$changelog" | grep -Po '^# \K.*')
    if [[ "$wdl_version" != "$changelog_version" ]]; then
      echo "MISMATCH: $(basename "$wdl") has version $wdl_version, changelog has $changelog_version"
    fi
  fi
done
```

**Common Causes of Smart Test Failure**:
- Merge conflicts resolved incorrectly
- Manual version updates applied to only one file
- Automated scripts updating WDL but not changelog (or vice versa)
- Cherry-picking commits that only update one file type
