# Whole Genome Germline Single Sample Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [WholeGenomeGermlineSingleSample_v2.2.0](https://github.com/broadinstitute/warp/releases) | January 6, 2021 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in WARP or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

The Whole Genome Germline Single Sample pipeline implements data pre-processing and initial variant calling according to the GATK Best Practices (June 2016) for germline SNP and Indel discovery in human whole-genome sequencing data. For a broad overview of the pipeline processes, read the GATK Best Practices documentation for [data pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) and for [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932).

The pipeline adheres to the Functional Equivalence pipeline specification ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)), a standard set of pipeline parameters to promote data interoperability across a multitude of global research projects and consortia. Read the [specification](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md) for full details or learn more about functionally equivalent pipelines in [this GATK blog](https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2018-02-09-Batch_effects_begone:_Introducing_the_Functional_Equivalence_data_processing_pipeline_spec.md).   

:::tip Want to try the Whole Genome Germline Single Sample pipeline?
You can test the pipeline in Terra! Go the [Whole-Genome-Analysis-Pipeline workspace](https://app.terra.bio/#workspaces/warp-pipelines/Whole-Genome-Analysis-Pipeline) which includes sample data and workflows for preprocessing and initial variant calling, sample map generation, and joint genotyping.
:::

## Set-up

### Workflow Installation and Requirements

The [Whole Genome Germline Single Sample workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl) is written in the Workflow Description Language [WDL](https://openwdl.org/) and can be downloaded by cloning the [warp repository](https://github.com/broadinstitute/warp/tree/master) in GitHub. The workflow can be deployed using [Cromwell](https://github.com/broadinstitute/cromwell), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. For the latest workflow version and release notes, please see the Whole Genome Germline Single Sample [changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.changelog.md).

### Software Version Requirements

* [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0) and GATK 3.5 for Whole Genome Variant Calling
* Picard 2.23.8
* Samtools 1.11
* Python 3.0
* Cromwell version support
    * Successfully tested on v52
    * Does not work on versions < v23 due to output syntax
* Papi version support
	* Successfully tested on Papi v2

### Input Requirements and Expectations

* Human whole-genome paired-end sequencing data in unmapped BAM (uBAM) format
* One or more read groups, one per uBAM file, all belonging to a single sample (SM)
* Input uBAM files must additionally comply with the following requirements:
    * All filenames have the same suffix (we use ".unmapped.bam")
    * Files pass validation by ValidateSamFile
    * Reads are in query-sorted order
    * All reads have an RG tag
* Reference genome must be Hg38 with ALT contigs

## Workflow Tasks and Tools

The Whole Genome Germline Single Sample [workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl) imports a series of tasks from the [tasks library](https://github.com/broadinstitute/warp/tree/master/tasks/broad) and a DNASeq struct ([DNASeqStructs.wdl](https://github.com/broadinstitute/warp/blob/master/structs/dna_seq/DNASeqStructs.wdl)) containing reference files from the [structs library](https://github.com/broadinstitute/warp/tree/master/structs).

Learn more about the software tools implemented in these tasks by reading the GATK [data pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) and [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) documentation.

:::tip Want to use the Whole Genome Germline Single Sample workflow in your publication?
Check out the workflow [Methods](./wgs.methods.md) to get started!
:::  


## Workflow Outputs

* CRAM, CRAM index, and CRAM MD5
* GVCF and its GVCF index
* BQSR report
* Summary metrics; to read more about any particular metric, you can search the metric using the [GATK documentation search](https://gatk.broadinstitute.org/hc/en-us/categories/360002302312)

### Base quality scores
The final CRAM files have base quality scores binned according to the [Functional Equivalence specification](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md#base-quality-score-binning-scheme) ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)).

| Original Score | Score after BQSR recalibration |
| --- | --- |
| 1-6 | unchanged |
| 7-12 | 10 | 
| 13-22 | 20 |
| 22-infinity | 30 |


## Important Notes

* Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
* For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
* Please visit the [GATK Technical Documentation](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591) site for further documentation on our workflows and tools.
* You can access relevant reference and resource bundles in the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811).

## Contact Us

This material is provided by the Data Sciences Platform group at the Broad Institute. Please direct any questions or concerns to one of our forum sites: [GATK](https://gatk.broadinstitute.org/hc/en-us/community/topics) or [Terra](https://support.terra.bio/hc/en-us/community/topics/360000500432).

## Licensing

Copyright Broad Institute, 2020 | BSD-3

This workflow is released under the WDL open source code license (BSD-3) (full license text at https://github.com/broadinstitute/warp/blob/master/LICENSE). However, please note that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.

Programs include:

- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml#13)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/terms/)
