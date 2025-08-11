# WARP AI Coding Guidelines

## Architecture Overview

WARP is a collection of cloud-optimized WDL (Workflow Description Language) pipelines for biological data processing. The repository follows a structured pattern with a **flattened directory structure**:

- **`pipelines/wdl/`**: All workflow definitions in a unified directory structure
- **`tasks/wdl/`**: All reusable WDL tasks in a unified directory structure  
- **`structs/`**: WDL struct definitions for type safety
- **`verification/`**: Test workflows that validate pipeline outputs
- **`scripts/`**: Build and validation automation

## WDL Development Patterns

### Import Structure
All WDL files use relative imports with this hierarchy:
```wdl
import "../../../tasks/wdl/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../structs/dna_seq/DNASeqStructs.wdl"
```

### Flattened Directory Structure
As of the latest reorganization, the repository uses a flattened structure:
- **Old structure**: `pipelines/broad/` and `pipelines/skylab/` → **New**: `pipelines/wdl/`
- **Old structure**: `tasks/broad/` and `tasks/skylab/` → **New**: `tasks/wdl/`

This maintains the same directory depth while consolidating team-specific directories.

### Pipeline Versioning
Every pipeline has a semantic version string that MUST be incremented for any changes:
```wdl
workflow MyPipeline {
  String pipeline_version = "1.2.4"  # Increment for changes
}
```

### Changelog Requirements
Each pipeline has a `<PipelineName>.changelog.md` file. For every version bump:
1. Add entry at top with format: `# 1.2.5\n2025-MM-DD (Date of Last Commit)\n\n* Description of changes`
2. Use past tense, active voice: "Updated docker image" not "Docker image was updated"

## Critical Workflows

### Validation
Always validate WDL syntax before committing:
```bash
java -jar ~/womtool-90.jar validate path/to/pipeline.wdl
```

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

## Team History and Pipeline Organization

- **Historical Teams**: Content originally organized by Broad and Skylab teams
  - Broad team: DNA-seq, arrays, reprocessing 
  - Skylab team: Single-cell RNA-seq, ATAC-seq, multiome
- **Current Organization**: Unified under `pipelines/wdl/` and `tasks/wdl/` directories
- Docker images maintained separately in [warp-tools](https://github.com/broadinstitute/warp-tools)

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
- Tests are triggered when files change in relevant `pipelines/wdl/` or `tasks/wdl/` paths
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

### Directory Structure
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
