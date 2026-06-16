---
title: Challenging Medically Relevant Genes Variant Calling Workflows (WDL)
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">

This section documents the All of Us CMRG workflows for remapping/calling in difficult regions and reblocking gVCFs for downstream joint analysis.

## Quick Summary

* **Purpose:** Improve variant calling in challenging medically relevant genes and prepare outputs for downstream aggregation.
* **Primary Outputs:** Corrected-region VCF/gVCF files and reblocked gVCFs.

## Workflow Overviews

| Step | Workflow | Description | WDL |
| :--: | --- | --- | :--: |
| 1 | [FixItFelix and Variant Calling](./fixitfelix_and_variant_call) | Subsets CMRG-relevant reads, remaps to masked GRCh38 using FixItFelix, and runs HaplotypeCaller. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.wdl) |
| 2 | [GVS AoU Reblock gVCF](./gvs_aou_reblock_gvcf) | Reblocks per-sample gVCFs and optionally copies outputs into AoU site-specific research buckets. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/GvsAoUReblockGvcf.wdl) |

## External Downstream Joint Calling

Joint calling can be performed with external GVS workflows, such as:

- [`GvsJointVariantCalling.wdl`](https://github.com/broadinstitute/gatk/blob/gvs_0.6.4/scripts/variantstore/wdl/GvsJointVariantCalling.wdl)

## Suggested Usage Pattern

1. Run FixItFelix remapping + variant calling (`FixItFelixAndVariantCall`).
2. Reblock gVCFs (`GvsAoUReblockGvcf`) as needed for storage/joint analysis.
3. Ingest and joint-call using external GVS pipelines.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>

