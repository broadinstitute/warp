---
sidebar_position: 4
---

# Version and Release Pipelines

WARP Pipelines are versioned semantically, allowing you to determine how and when your data was created (provenance). This promotes compatibility across datasets and ensures that analyses can be reproduced by the global scientific community. Semantic version numbers (written as major.minor.patch) are human readable and give immediate insight into the compatibility of pipeline outputs (see the [Versioning Guidelines](#versioning-guidelines) below).

Versions of each pipeline are packaged into releases and published on GitHub (see the WARP [releases page](https://github.com/broadinstitute/warp/releases)). A published release of a pipeline version in GitHub has passed scientific tests (read more in [TestingPipelines](./TestingPipelines.md)) and is available to be used in production systems.

:::tip To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/develop/wreleaser).
:::

## Versioning Requirements

All versioned pipelines must have:
* A release version number in the main workflow WDL. The version number is in the form `major.minor.patch` and is included as a field in the standardized format `String pipeline_version = "1.0.0"`.
* A cumulative changelog file containing the changes introduced in the PR and all previous changes in reverse chronological order. Changelog entries are to be formatted according to the [changelog style guide](../contribution/contribute_to_warp/changelog_style.md) and will include the version number, date of last commit, and a bulleted list of changes since the last release.

## Versioning Guidelines

In WARP, a pipeline requires a version change when any change is made to the pipeline’s main WDL workflow or any of the WDL workflow dependencies. If a change requires a new pipeline version number, the changes and new version number are demarcated in the pipeline’s changelog.

Pipeline version numbers are updated based the following  guidelines:
#### Major changes
* Any qualitative change to the pipeline’s scientific outputs. If you use the pipeline’s data output, this change indicates a possible need to reprocess data analyzed with a previous release version.
* Any breaking changes to the pipeline, including input/output refactors, renaming of the pipeline, and changes to input/output formats.

#### Minor changes
* Addition of new outputs that don’t impact previous outputs; for example, adding a new md5 checksum output or outputting new QC metrics.
* Changes to the pipeline that do not qualitatively impact the scientific outputs, but may produce slightly different outputs (no data reprocessing needed).

#### Patch (micro) changes
* Memory changes, internal refactor or variable name changes, speed or cost optimizations, comments, metadata.
* Addition of optional inputs.

When pipelines are promoted to the master branch, a script packages the pipeline for release on GitHub. A release contains three components:
1. A release name comprising the pipeline name and version number listed in the changelog (i.e. SmartSeq2SingleSample_v5.0.0)
2. Release notes comprising the corresponding version changelog entry
3. Artifacts including the main workflow WDL, a zip of all workflow dependencies, and when applicable, an options file

Upon release, the pipeline is automatically pushed to Dockstore based on the WARP [Dockstore configuration](https://github.com/broadinstitute/warp/blob/develop/.dockstore.yml).

To support early integration testing of our pipelines, we also maintain a floating “dev” pre-release for each pipeline named “PipelineName_develop”.
