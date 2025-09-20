---
sidebar_position: 1
slug: /Pipelines/Exome_Germline_Single_Sample_Pipeline/README
---

# Exome Germline Single Sample Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [ExomeGermlineSingleSample_v3.1.19](https://github.com/broadinstitute/warp/releases?q=ExomeGermlineSingleSample_v3.0.0&expanded=true) | March, 2024 | Elizabeth Kiernan | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues). |


The Exome Germline Single Sample pipeline implements data pre-processing and initial variant calling according to the GATK Best Practices for germline SNP and Indel discovery in human exome sequencing data.

For a broad overview of the pipeline processes, read the GATK Best Practices documentation for [data pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) and for [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932).

:::tip Want to try the Exome Germline Single Sample pipeline?
You can test the pipeline in Terra! Go the [Exome-Analysis-Pipeline workspace](https://app.terra.bio/#workspaces/warp-pipelines/Exome-Analysis-Pipeline) which includes sample data and workflows for preprocessing and initial variant calling, sample map generation, and joint genotyping.
:::


## Set-up

### Workflow Installation and Requirements

The Exome Germline Single Sample workflow is written in the Workflow Description Language WDL and can be downloaded by cloning the [warp repository](https://github.com/broadinstitute/warp/tree/master) in GitHub. The workflow can be deployed using [Cromwell](https://github.com/broadinstitute/cromwell), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. For the latest workflow version and release notes, please see the Exome Germline Single Sample [changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.changelog.md).

### Software Version Requirements

* [GATK 4.5.0.0](https://github.com/broadinstitute/gatk/releases/tag/4.5.0.0)
* Picard 2.26.10
* Samtools 1.11
* Python 3.0
* Cromwell version support
    * Successfully tested on v52
    * Does not work on versions < v23 due to output syntax
* Papi version support
	* Successfully tested on Papi v2

### Input Requirements and Expectations

- Human exome sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
    * Filenames all have the same suffix (we use ".unmapped.bam")
    * Files must pass validation by ValidateSamFile
    * Reads are provided in query-sorted order
    * All reads must have an RG tag
- Reblocked GVCF output names must end in ".rb.g.vcf.gz"
- Reference genome must be Hg38 with ALT contigs
- Unique exome calling, target, and bait [.interval_list](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852) obtained from the sequencing provider. Generally the calling, target, and bait files will not be the same.


## Workflow Tasks and Tools

The Exome Germline Single Sample [workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl) imports a series of tasks from the WARP [tasks library](https://github.com/broadinstitute/warp/tree/master/tasks/wdl) and a DNASeq struct ([DNASeqStructs.wdl](https://github.com/broadinstitute/warp/blob/master/structs/dna_seq/DNASeqStructs.wdl)) containing reference files from the [structs library](https://github.com/broadinstitute/warp/tree/master/structs/dna_seq).

You can read more about the software tools implemented in these tasks by reading the GATK [data pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) and [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) documentation.

:::tip Want to use the Exome Germline Single Sample workflow in your publication?
Check out the workflow [Methods](./exome.methods.md) to get started!
:::

## Workflow Outputs

- CRAM, CRAM index, and CRAM MD5
- [Reblocked](https://gatk.broadinstitute.org/hc/en-us/articles/360037593171) GVCF and its GVCF index (read more in the [Reblocking](#reblocking) section below)
- BQSR report
- Summary metrics; to read more about any particular metric, you can search the metric using the [GATK documentation search](https://gatk.broadinstitute.org/hc/en-us/categories/360002302312)

### Reblocking
Reblocking is a process that compresses a HaplotypeCaller GVCF by merging homRef blocks according to new genotype quality (GQ) bands and facilitates joint genotyping by removing alt alleles that do not appear in the called genotype. 

As of November 2021, reblocking is a default task in the Exome pipeline. To skip reblocking, add the following to the workflow's input configuration file (JSON):

```WDL
"ExomeGermlineSingleSample.BamToGvcf.skip_reblocking": true
```

The [Reblocking task](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/GermlineVariantDiscovery.wdl) uses the GATK ReblockGVCF tool with the arguments:

```WDL
-do-qual-approx -floor-blocks -GQB 20 -GQB 30 -GQB 40 
```
The following summarizes how reblocking affects the Exome GVCF and downstream tools compared to the GVCF produced with the default HaplotypeCaller GQ bands:


1. PLs are omitted for homozygous reference sites to save space– GQs are output for genotypes, PLs can be approximated as [0, GQ, 2\*GQ].

2. GQ resolution for homozygous reference genotypes is reduced (i.e. homRef GQs will be underconfident) which may affect analyses like de novo calling where well-calibrated reference genotype qualities are important.

3. Alleles that aren’t called in the sample genotype are dropped. Each variant should have no more than two non-symbolic alt alleles, with the majority having just one plus <NON_REF>.

4. New annotations enable merging data for filtering without using genotypes. For example:
    * RAW_GT_COUNT(S) for doing ExcessHet calculation from a sites-only file.
    * QUALapprox and/or AS_QUALapprox for doing QUAL approximation/filling. 
    * QUAL VCF field from a combined sites-only field.
    * VarDP and/or AS_VarDP used to calculate QualByDepth/QD annotation for VQSR.

5. The MIN_DP has been removed.

6. Reblocked GVCFs have the following cost/scale improvements:
    * A reduced storage footprint compared with HaplotypeCaller GVCF output.
    * Fewer VariantContexts (i.e. lines) per VCF which speeds up GenomicsDB/Hail import.
    * Fewer alternate alleles which reduce memory requirements for merging.

Additionally, the 4 GQ band schema has specific improvements compared with the 7-band schema:
1. It does not drop GQ0s; reblocked GVCFs should cover all the positions that the input GVCF covers.
2. It has no overlaps; the only overlapping positions should be two variants (i.e. deletions) on separate haplotypes.
3. No more no-calls; all genotypes should be called. Positions with no data will be homRef with GQ0.

Read more about the reblocked GVCFs in the [WARP Blog](https://broadinstitute.github.io/warp/blog/Nov21_ReblockedGVCF).


### Base quality scores
The final CRAM files have base quality scores binned according to the [Functional Equivalence specification](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md#base-quality-score-binning-scheme) ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)).

| Original Score | Score after BQSR recalibration |
| --- | --- |
| 1-6 | unchanged |
| 7-12 | 10 |
| 13-22 | 20 |
| 22-infinity | 30 |

## Important Notes

- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Please visit the [GATK Technical Documentation](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591) site for further documentation on our workflows and tools.
- You can access relevant reference and resource bundles in the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811).

## Citing the Exome Germline Single Sample Pipeline

If you use the Exome Germline Single Sample Pipeline in your research, please consider citing our preprint:

Degatano, K., Awdeh, A., Cox III, R.S., Dingman, W., Grant, G., Khajouei, F., Kiernan, E., Konwar, K., Mathews, K.L., Palis, K., et al. Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis. Bioinformatics 2025; btaf494. https://doi.org/10.1093/bioinformatics/btaf494

## Contact Us

This material is provided by the Data Science Platform group at the Broad Institute. Please direct any questions or concerns to one of our forum sites : [GATK](https://gatk.broadinstitute.org/hc/en-us/community/topics) or [Terra](https://support.terra.bio/hc/en-us/community/topics/360000500432).

## Licensing

Copyright Broad Institute, 2023 | BSD-3

The workflow script is released under the **WDL open source code license (BSD-3)** (full license text at https://github.com/broadinstitute/warp/blob/master/LICENSE). However, please note that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.

- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml#13)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/terms/)

