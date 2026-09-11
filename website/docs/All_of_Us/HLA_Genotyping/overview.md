---
title: HLA Genotyping Workflows (WDL)
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">
This section documents the All of Us HLA workflows for per-sample HLA calling and cohort-level aggregation.

## Quick Summary

* **Purpose:** Perform HLA genotyping from aligned sequencing data and aggregate per-sample consensus calls.
* **Primary Outputs:** Per-sample consensus HLA calls and combined cohort-level tables.

## Workflow Overviews

| Step | Workflow | Description | WDL |
| :--: | --- | --- | :--: |
| 1 | [HLA Consensus Genotyping](./hla_consensus_genotyping) | Extracts HLA reads, runs HLA-HD with conditional Polysolver/OptiType fallback, and emits consensus calls. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping.wdl) |
| 2 | [HLA Make Table](./make_table) | Combines per-sample consensus calls into a cohort-level result table. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/MakeTable.wdl) |

## Suggested Usage Pattern

1. Run `HLAGenotyping` per sample.
2. Collect consensus outputs and run `MakeTable` for cohort aggregation.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
