# WARP AI Coding Guidelines

## Architecture Overview

WARP is a collection of cloud-optimized WDL (Workflow Description Language) pipelines for biological data processing. The repository follows a structured pattern:

- **`pipelines/wdl/`**: Main workflow definitions organized by domain (broad vs skylab teams)
- **`tasks/`**: Reusable WDL tasks imported by pipelines 
- **`structs/`**: WDL struct definitions for type safety
- **`verification/`**: Test workflows that validate pipeline outputs
- **`scripts/`**: Build and validation automation

## WDL Development Patterns

### Import Structure
All WDL files use relative imports with this hierarchy:
```wdl
import "../../../tasks/broad/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../structs/dna_seq/DNASeqStructs.wdl"
```

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

## Team-Specific Patterns

- **Broad team**: DNA-seq, arrays, reprocessing (pipelines/wdl/broad/)
- **Skylab team**: Single-cell RNA-seq, ATAC-seq, multiome (pipelines/wdl/skylab/)
- Docker images maintained separately in [warp-tools](https://github.com/broadinstitute/warp-tools)

## Key Files
- `pipeline_versions.txt`: Master version tracking
- `changelog_style.md`: Detailed changelog formatting rules  
- `wreleaser/`: CLI tool for querying releases
- `scripts/common.sh`: Shared build functions
