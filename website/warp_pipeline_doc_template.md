---
sidebar_position: 1
slug: /Pipelines/<PipelineName>/README
---

# <Pipeline Name> Overview

<!-- 
  METADATA TABLE
  - Pipeline Version: current WDL version string (e.g., "1.0.0")
  - Date Updated: human-readable date (e.g., "March, 2025")
  - Documentation Author: full name or team name
  - Questions or Feedback: link to WARP GitHub issues
-->

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [v1.0.0](https://github.com/broadinstitute/warp/releases) | March, 2025 | [WARP Pipelines](mailto:warp@broadinstitute.org) | [File an issue](https://github.com/broadinstitute/warp/issues) |

<!-- 
  OPTIONAL: Pipeline diagram image
  Omit if no diagram exists.
  Example:
    ![<Pipeline Name> diagram](./diagram.png)
-->

## Introduction to the <Pipeline Name> workflow

<!-- 
  2–4 paragraph narrative describing:
  1. What the pipeline does at a high level and what data type it processes.
  2. The main steps (e.g., alignment, duplicate marking, QC, output generation).
  3. Any major variants or modes the pipeline supports (e.g., DRAGEN mode, sc vs. sn).
  4. Links to any relevant count matrix or metrics overviews if applicable.
-->

<Pipeline Name> is a [cloud-optimized](https://github.com/broadinstitute/warp) WDL pipeline that processes [data type] data. The workflow [brief description of main processing steps].

[Additional paragraph about key features, tools used, or output formats.]

[Optional: Paragraph about any pipeline modes, consortium use cases, or special configurations.]

## Quickstart table

<!-- 
  3-column summary table for users to quickly assess pipeline suitability.
  Common rows: assay type, cell/sample scope, supported species, input data type,
  output file formats, workflow language, technology support.
  OMIT this section for array-based or non-single-cell pipelines (e.g., Illumina Genotyping Arrays,
  Joint Genotyping) that do not follow the single-cell pattern.
-->

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Assay type | [e.g., Single-nucleus RNA-seq] | [Link to assay info] |
| Overall workflow | [e.g., Quality control, alignment, quantification] | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic reference sequence | GRCh38 / mm10 | [GENCODE](https://www.gencodegenes.org/) |
| Aligner | [e.g., STARsolo] | [Link] |
| Data input file format | FASTQ | |
| Data output file format | [e.g., BAM, h5ad, CSV] | |

## Set-up

### <Pipeline Name> installation and requirements

The <Pipeline Name> workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [<Pipeline Name> changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/<pipeline_dir>/<PipelineName>.changelog.md).

<!--
  Changelog link pattern:
  - Most pipelines: pipelines/wdl/<dir>/<PipelineName>.changelog.md
  - Some use underscore: pipelines/wdl/<dir>/<PipelineName>_changelog.md
  - wreleaser mention (include for all pipelines):
-->

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system. To learn more about Cromwell, see the [Cromwell documentation](https://cromwell.readthedocs.io/en/stable/).

<!--
  OPTIONAL: Terra workspace tip block.
  Include only if a Terra workspace exists for this pipeline.

:::tip
You can also try this pipeline in Terra! See the [<Pipeline Name> Terra workspace](https://app.terra.bio/#workspaces/...) for a preconfigured environment with sample data.
:::
-->

## Inputs

<!--
  Input description guidance:
  - Simple pipelines (≤20 inputs): use a single flat table.
  - Complex pipelines (many inputs): use typed sub-tables grouped by category
    (e.g., Sample inputs, Reference inputs, Optional inputs, Contamination inputs).
  - Column options: "Input variable name | Description | Type" (3-col)
                    "Input variable name | Type | Description" (3-col, common)
                    "Input variable name | Type | Description | Default value" (4-col for optional inputs)
  - Struct inputs: note the struct name and link to the struct WDL definition.
-->

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `input_variable_1` | Description of the first input. | String |
| `input_variable_2` | Description of the second input. | File |
| `input_variable_3` | *(Optional)* Description with default value. Default: `false`. | Boolean |

<!--
  If the pipeline uses WDL structs, add a sentence like:
  "Some inputs are defined using the `<StructName>` struct (see [`<StructName>.wdl`](link))."
-->

## <Pipeline Name> tasks and tools

<!--
  Opening paragraph pattern (use verbatim or adapt):
  "The <Pipeline Name> pipeline calls a series of tasks defined in the [<PipelineName>.wdl](link) and [tasks WDL](link)."
  Then list the overall steps as a numbered list, followed by the task table.
-->

The <Pipeline Name> pipeline calls a series of tasks to [brief description]. The following sections describe each step:

1. [Step 1 description](#section-anchor)
2. [Step 2 description](#section-anchor)
3. [Step N description](#section-anchor)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `String docker =`.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [TaskName](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/<dir>/<Pipeline>.wdl) | Tool name | Software name | Description of what this task does. |

### 1. Step 1 name

<!-- Describe the step's purpose, inputs consumed, tools used, and outputs produced. -->

### 2. Step 2 name

<!-- Repeat for each major step. -->

## Outputs

<!--
  Output table: 2 or 3 columns.
  Common patterns:
    "Output variable name | Filename, if applicable | Output format and description"  (3-col, single-cell)
    "Output variable name | Description | Type"  (3-col, genomics)
    "Output variable name | Description"  (2-col, simpler pipelines)
-->

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `output_variable_1` | `<input_id>.bam` | BAM file containing [description]. |
| `output_variable_2` | `<input_id>.h5ad` | h5ad (Anndata) containing [description]. |

## Versioning and testing

<!--
  Standard versioning blurb — adapt pipeline name and changelog path.
  Use "Versioning and testing" for single-cell pipelines; "Versioning" is also common.
-->

All <Pipeline Name> pipeline releases are documented in the [<Pipeline Name> changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/<pipeline_dir>/<PipelineName>.changelog.md) and tested using [plumbing and scientific test data](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/<pipeline_dir>/test_inputs). To learn more about WARP pipeline testing, see [Testing Pipelines](https://broadinstitute.github.io/warp/docs/About_WARP/TestingPipelines).

## Citing the <Pipeline Name> Pipeline

<!--
  Standard WARP citation block. Use verbatim. Optionally precede with a SciCrunch RRID block
  if this pipeline has a registered resource identifier.
-->

<!--
  OPTIONAL: SciCrunch RRID block (include only if an RRID has been registered):

If you use the <Pipeline Name> Pipeline in your research, please identify the pipeline in your methods section using the [<Pipeline Name> SciCrunch resource identifier](https://scicrunch.org/resources).

* Ex: *<Pipeline Name> Pipeline (RRID:SCR_XXXXXX)*
-->

When citing WARP, please use the following:

Kylee Degatano, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, Kishori Konwar, Kaylee L Mathews, Kevin Palis, Nikelle Petrillo, Geraldine Van der Auwera, Chengchen (Rex) Wang, Jessica Way. "Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis." _Bioinformatics_, 2025; [https://doi.org/10.1093/bioinformatics/btaf494](https://doi.org/10.1093/bioinformatics/btaf494).

## Acknowledgements

<!--
  OPTIONAL section. Include when there are specific external collaborators or consortia to credit.
  Omit for pipelines without notable external contributions.
  Standard opening: "We are immensely grateful to the members of [consortium/initiative]..."
-->

We are immensely grateful to the members of [Consortium or Initiative] for their invaluable contributions to this pipeline.

## Important notes

<!--
  OPTIONAL section used in genomics pipelines (WGS, Exome, etc.) for caveats,
  runtime optimization notes, and licensing reminders.
  Common entries:
  - "Runtime parameters are optimized for Broad's Google Cloud Platform implementation."
  - Links to GATK documentation.
  - DRAGEN-specific notes.
  Omit for single-cell pipelines where this section does not typically appear.
-->

- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.

## Licensing

<!--
  OPTIONAL but recommended for genomics pipelines. Omit for single-cell pipelines
  where this section is not standard.
-->

Copyright Broad Institute, 2025 | BSD-3

The workflow script is released under the **WDL open source code license (BSD-3)** (full license text at https://github.com/broadinstitute/warp/blob/master/LICENSE). However, please note that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.

## Feedback

<!--
  Standard feedback section — use verbatim.
  Some pipelines use "Contact us" instead of "Feedback" — either is acceptable.
-->

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

<!--
  OPTIONAL: FAQs section — include only for high-traffic pipelines with known user questions
  (e.g., Optimus, WGS). Use Docusaurus :::note blocks:

## FAQs

:::note Question Can I run <Pipeline Name> in Terra?

Yes! [Answer here.]
:::
-->
