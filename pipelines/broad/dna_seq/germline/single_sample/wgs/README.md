| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [WholeGenomeGermlineSingleSample_v2.0](https://github.com/broadinstitute/warp/releases) | October 02, 2020 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in WARP or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

# Introduction to the Whole Genome Germline Single Sample Pipeline

The Whole Genome Germline Single Sample pipeline implements data pre-processing and initial variant calling according to the GATK Best Practices (June 2016) for germline SNP and Indel discovery in human whole-genome sequencing data. For a broad overview of the pipeline processes, read the GATK Best Practices documentation for [data pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) and for [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932).

# Set-up

## Workflow Installation and Requirements

The [Whole Genome Germline Single Sample workflow](WholeGenomeGermlineSingleSample.wdl) is written in the Workflow Description Language [WDL](https://openwdl.org/) and can be downloaded by cloning the GitHub repository [WARP](https://github.com/broadinstitute/warp/). The workflow can be deployed using [Cromwell](https://github.com/broadinstitute/cromwell), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. For the latest workflow version and release notes, please see the Whole Genome Germline Single Sample [changelog](WholeGenomeGermlineSingleSample.changelog.md).

## Software Version Requirements
* GATK 3.5 and GATK 4.beta.5
* Picard 2.20.0-SNAPSHOT
* Samtools 1.3.1
* Python 2.7 and 3.0
* Cromwell version support 
    * Tested on Cromwell 52
    * Does not work on versions < v23 due to output syntax


## Input Requirements and Expectations

* Human whole-genome paired-end sequencing data in unmapped BAM (uBAM) format
* One or more read groups, one per uBAM file, all belonging to a single sample (SM)
* Input uBAM files must additionally comply with the following requirements:
    * All filenames have the same suffix (we use ".unmapped.bam")
    * Files pass validation by ValidateSamFile
    * Reads are in query-sorted order
    * All reads have an RG tag
* Reference genome must be Hg38 with ALT contigs

# Workflow Tasks and Tools

The Whole Genome Germline Single Sample [workflow](WholeGenomeGermlineSingleSample.wdl) imports a series of tasks from the [tasks library](../../../../../../tasks/broad/) and a DNASeq struct ([DNASeqStructs.wdl](../../../../../../structs/dna_seq/DNASeqStructs.wdl)) containing reference files from the [structs library](../../../../../../structs/).

You can read more about the software tools implemented in these tasks by reading the GATK [data pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) and [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) documentation.

# Workflow Outputs
* CRAM, CRAM index, and CRAM MD5
* GVCF and its GVCF index
* BQSR report
* Summary metrics; to read more about any particular metric, you can search the metric using the [GATK documentation search](https://gatk.broadinstitute.org/hc/en-us/categories/360002302312)


# Important Notes
* The accompanying JSON is a generic, ready to use, example template for the workflow. It is the user’s responsibility to correctly set the reference and resource variables for their own particular test case using the [GATK Tool and Tutorial Documentations](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591).
* By default, HaplotypeCaller will perform variant calling using GATK 3.5, which is what is used in Broad Production. To use GATK4, specify `use_gatk3_haplotype_caller=false` in the inputs.json.
* Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
* For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
* Please visit the [GATK Technical Documentation](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591) site for further documentation on our workflows and tools.
* You can access relevant reference and resource bundles in the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811).

# Contact Us 
This material is provided by the Data Sciences Platform group at the Broad Institute. Please direct any questions or concerns to one of our forum sites: [GATK](https://gatk.broadinstitute.org/hc/en-us/community/topics) or [Terra](https://support.terra.bio/hc/en-us/community/topics/360000500432).

# Licensing
Copyright Broad Institute, 2020 | BSD-3

This workflow is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). However, please note that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. Programs include:
- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml#13)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/terms/)

