# WARP Agents and Tools

This document provides an overview of agents, tools, and utilities available in the WARP repository, organized by function and location.

## Table of Contents

- [Build and Release Tools](#build-and-release-tools)
- [WDL Processing Tools](#wdl-processing-tools)
- [Data Processing Scripts](#data-processing-scripts)
- [API Clients](#api-clients)
- [Code Generation and Analysis](#code-generation-and-analysis)
- [Environment and Terminal Execution Guide](#environment-and-terminal-execution-guide)
- [How to Add New Tools](#how-to-add-new-tools)

---

## Build and Release Tools

### Wreleaser

**Location:** [`wreleaser/`](wreleaser/)

A command-line tool for managing and querying WARP pipeline releases. Use Wreleaser to:
- Search and discover available pipeline releases
- Retrieve release information and metadata
- Automate release management workflows

**Related Files:**
- Main binary and source code in `wreleaser/`

### Build and Release Scripts

**Location:** [`scripts/`](scripts/)

Automated scripts for building, validating, and releasing pipelines:

- [build_pipeline_release.sh](scripts/build_pipeline_release.sh) - Build pipeline releases with version management
- [build_workflow_release.sh](scripts/build_workflow_release.sh) - Build individual workflow releases
- [release_pipeline_to_github.sh](scripts/release_pipeline_to_github.sh) - Release pipelines to GitHub
- [production_release.sh](scripts/production_release.sh) - Production release automation
- [validate_wdls.sh](scripts/validate_wdls.sh) - Validate WDL syntax and correctness
- [validate_release.sh](scripts/validate_release.sh) - Validate release artifacts

**Usage:**
```bash
cd scripts/
./build_pipeline_release.sh -w pipeline.wdl -v version -o output_dir
```

---

## WDL Processing Tools

### WDL to Bash Converter

**Location:** [`scripts/wdl_to_bash/`](scripts/wdl_to_bash/)

Tools for extracting and converting WDL task definitions to executable bash scripts, including a test mode for bash error detection.

**Files:**
- [wdl_to_bash_extractor.py](scripts/wdl_to_bash/wdl_to_bash_extractor.py) - Main converter utility with test mode
- [test_wdl_to_bash_extractor.py](scripts/wdl_to_bash/test_wdl_to_bash_extractor.py) - Test suite
- [example_extraction_workflow.sh](scripts/wdl_to_bash/example_extraction_workflow.sh) - Example workflow

**Documentation:**
- [WDL_TO_BASH_README.md](scripts/wdl_to_bash/WDL_TO_BASH_README.md) - Getting started guide
- [WDL_TO_BASH_USAGE_GUIDE.md](scripts/wdl_to_bash/WDL_TO_BASH_USAGE_GUIDE.md) - Detailed usage instructions
- [WDL_TO_BASH_ADVANCED.md](scripts/wdl_to_bash/WDL_TO_BASH_ADVANCED.md) - Advanced usage patterns
- [PACKAGE_SUMMARY.md](scripts/wdl_to_bash/PACKAGE_SUMMARY.md) - Package overview
- [QUICK_REFERENCE.md](scripts/wdl_to_bash/QUICK_REFERENCE.md) - Quick reference guide

**Usage - Extract Mode:**
```bash
python3 /workspaces/warp/scripts/wdl_to_bash/wdl_to_bash_extractor.py /workspaces/warp/pipeline.wdl -o /workspaces/warp/tmp/extracted_tasks
```

**Usage - Test Mode** (detect bash errors without extracting):
```bash
python3 /workspaces/warp/scripts/wdl_to_bash/wdl_to_bash_extractor.py /workspaces/warp/pipeline.wdl --test
```

Test mode identifies bash syntax errors in all task command blocks and generates a detailed error report without creating any output files. Use this for quality assurance and debugging before extraction.

---

## Data Processing Scripts

### Common Utilities

**Location:** [`scripts/`](scripts/)

Shared functions and utilities used by other scripts:

- [common.sh](scripts/common.sh) - Common functions and helper utilities

### Interval and Region Generation

**Location:** [`scripts/`](scripts/)

Utilities for generating genomic intervals and regions:

- [generate_hg38_exome_calling_intervals.v1.sh](scripts/generate_hg38_exome_calling_intervals.v1.sh) - Generate exome calling intervals (v1)
- [generate_hg38_exome_calling_intervals.v1.1.sh](scripts/generate_hg38_exome_calling_intervals.v1.1.sh) - Generate exome calling intervals (v1.1)
- [generate_bge_exome_calling_regions.sh](scripts/generate_bge_exome_calling_regions.sh) - Generate BGE exome calling regions
- [subsetWGSIntervals.README](scripts/subsetWGSIntervals.README) - Documentation for WGS interval subsetting

### Data Analysis Scripts

**Location:** [`scripts/`](scripts/)

Scripts for cost calculation and data analysis:

- [calculate_cost.py](scripts/calculate_cost.py) - Calculate pipeline execution costs
- [removeBadBlats.py](scripts/removeBadBlats.py) - Filter poor quality BLAT results

---

## API Clients

### Dockstore API

**Location:** [`scripts/dockstore_api/`](scripts/dockstore_api/)

Client library for interacting with the Dockstore registry API.

### FireCloud API

**Location:** [`scripts/firecloud_api/`](scripts/firecloud_api/)

Client library for interacting with the FireCloud (Terra) API.

### Pipeline Client

**Location:** [`scripts/`](scripts/)

- [pipeline_client.sh](scripts/pipeline_client.sh) - General pipeline client utilities

---

## Code Generation and Analysis

### WDL Analysis and Validation

**Location:** [`scripts/`](scripts/)

Tools for analyzing and validating WDL workflows:

- [get_changed_pipeline_worklow_test_args.sh](scripts/get_changed_pipeline_worklow_test_args.sh) - Get test arguments for changed pipelines
- [BuildAFComparisonTable.wdl](scripts/BuildAFComparisonTable.wdl) - Helper workflow for AF comparison
- [FilterAFComparisonTable.wdl](scripts/FilterAFComparisonTable.wdl) - Helper workflow for filtering AF data
- [RemoveBadSitesByID.wdl](scripts/RemoveBadSitesByID.wdl) - Helper workflow for site filtering

---

## Configuration and Documentation

- [.dockstore.yml](.dockstore.yml) - Dockstore registry configuration
- [.github/](github/) - GitHub workflows and automation
- [pipeline_versions.txt](pipeline_versions.txt) - Master version tracking for all pipelines
- [changelog_style.md](changelog_style.md) - Style guide for pipeline changelogs

---

## Environment and Terminal Execution Guide

### Dev Container Context

This repository operates within a Linux dev container (Ubuntu 22.04.5 LTS) with specific execution characteristics:

#### Path Resolution in Terminal Commands

**Issue Observed:** Relative paths in Python/shell command invocations can fail due to working directory context changes during command execution.

**Solution:** Always use **absolute paths** when invoking Python scripts and tools from the terminal:

```bash
# ❌ AVOID (may fail due to path resolution):
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py pipeline.wdl

# ✅ RECOMMENDED (reliable):
python3 /workspaces/warp/scripts/wdl_to_bash/wdl_to_bash_extractor.py /workspaces/warp/all_of_us/mitochondria/pipeline.wdl -o extracted_tasks
```

#### Best Practices for Agent Terminal Execution

When running terminal commands in this environment:

1. **Use absolute paths** for executables and input files
2. **Use absolute paths** for output directories (especially for Python scripts)
3. **Avoid changing working directories** unless necessary—prefer absolute paths instead
4. **Test with full paths first** if relative paths fail unexpectedly
5. **File location check**: Verify file existence before executing commands:
   ```bash
   ls -la /path/to/file  # Verify file exists
   python3 /absolute/path/to/script.py /absolute/path/to/input.wdl
   ```

#### WDL to Bash Extractor Execution

For the WDL to Bash extractor tool, always use absolute paths:

```bash
# Correct usage
/usr/bin/python3 /workspaces/warp/scripts/wdl_to_bash/wdl_to_bash_extractor.py \
  /workspaces/warp/all_of_us/mitochondria/mt_coverage_merge.wdl \
  -o extracted_mt_coverage -l INFO
```

**Why:** The extractor uses Python's file resolution internally, and relative paths can break when the terminal context changes during execution.

#### Temporary Files and Output Management

**Important:** Always store temporary files and intermediate outputs in the `tmp/` directory to prevent accidental commits to git.

```bash
# ✅ CORRECT - Uses tmp/ directory (excluded from git)
python3 /workspaces/warp/scripts/wdl_to_bash/wdl_to_bash_extractor.py \
  /workspaces/warp/all_of_us/mitochondria/pipeline.wdl \
  -o /workspaces/warp/tmp/extracted_tasks -l INFO

# ❌ AVOID - Creates files in repo that may be committed
python3 /workspaces/warp/scripts/wdl_to_bash/wdl_to_bash_extractor.py \
  /workspaces/warp/all_of_us/mitochondria/pipeline.wdl \
  -o extracted_tasks -l INFO
```

**Guidelines:**
- All test outputs, extracted files, and temporary results go in `/workspaces/warp/tmp/`
- The `tmp/` directory is listed in [.gitignore](.gitignore) and will not be committed
- Create subdirectories within `tmp/` for organization: `tmp/extracted_tasks/`, `tmp/test_results/`, etc.
- Clean up large temporary files after use to avoid filling disk space
- If generating test data or verification outputs, use `verification/` directory instead (for committed test fixtures)

---

## How to Add New Tools

When adding new tools or agents to WARP:

1. **Script-based tools** go in [`scripts/`](scripts/) or a subdirectory if related (e.g., `scripts/wdl_to_bash/`)
2. **Client libraries** go in subdirectories like [`scripts/firecloud_api/`](scripts/firecloud_api/)
3. **Documentation** accompanies the tools in their respective directories
4. **Top-level documentation** (like this file) provides navigation and overview

---

## Further Reading

- [WARP WDL Style Guide](WARP_WDL_Style_Guide.md)
- [README.md](README.md) - Main project documentation
- [WARP Documentation Site](https://broadinstitute.github.io/warp/)
