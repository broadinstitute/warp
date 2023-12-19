---
sidebar_position: 1
slug: /Pipelines/BuildIndices_Pipeline/README
---

# BuildIndices Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [BuildIndices_v3.0.0](https://github.com/broadinstitute/warp/releases) | December, 2023 | Kaylee Mathews | Please file GitHub issues in warp or contact [documentation authors](mailto:warp-pipelines-help@broadinstitute.org) |

![BuildIndices_diagram](./buildindices_diagram.png)


## Introduction to the BuildIndices workflow

The [BuildIndices workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/build_indices/BuildIndices.wdl) is an open-source, cloud-optimized pipeline developed in collaboration with the [BRAIN Initiative Cell Census Network](https://biccn.org/) (BICCN) and the BRAIN Initiative Cell Atlas Network (BICAN). 

Overall, the workflow filters GTF files for selected gene biotypes, calculates chromosome sizes, and builds reference bundles with required files for [STAR](https://github.com/alexdobin/STAR) and [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) aligners.

## Quickstart table
The following table provides a quick glance at the BuildIndices pipeline features:

| Pipeline features | Description | Source |
| --- | --- | --- |
| Overall workflow | Reference bundle creation for STAR and bwa-mem2 aligners | Code available on [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/build_indices/BuildIndices.wdl) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence, M32 (GRCm39) mouse genome primary sequence, and release 103 (GCF_003339765.1) macaque genome primary sequence  | GENCODE [human reference files](https://www.gencodegenes.org/human/release_43.html), GENCODE [mouse reference files](https://www.gencodegenes.org/mouse/release_M32.html), and NCBI [macaque reference files](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003339765.1/) |
| Gene annotation reference (GTF) | Reference containing gene annotations | GENCODE [human GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz), GENCODE [mouse GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz), and NCBI [macaque GTF](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003339765.1/) |
| Reference builders | STAR, bwa-mem2 | [Dobin et al. 2013](https://pubmed.ncbi.nlm.nih.gov/23104886/), [Vasimuddin et al. 2019](https://ieeexplore.ieee.org/document/8820962) |
| Data input file format | File format in which reference files are provided | FASTA, GTF, TSV |
| Data output file format | File formats in which BuildIndices output is provided | GTF, TAR, TXT |

## Set-up

### BuildIndices installation

To download the latest BuildIndices release, see the release tags prefixed with "BuildIndices" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All BuildIndices pipeline releases are documented in the [BuildIndices changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/build_indices/BuildIndices.changelog.md). 

To search releases of this and other pipelines, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

If you’re running a BuildIndices workflow version prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the folder `website/docs/Pipelines/BuildIndices_Pipeline`).

The BuildIndices pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform.

### Inputs

The BuildIndices workflow inputs are specified in JSON configuration files. Configuration files for [macaque](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/build_indices/Macaque.json) and [mouse](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/build_indices/Mouse.json) references can be found in the WARP repository.

#### Input descriptions

| Parameter name | Description | Type |
| --- | --- | --- |
| genome_source | Describes the source of the reference genome listed in the GTF file; used to name output files; can be set to “NCBI” or “GENCODE”. | String |
| gtf_annotation_version | Version or release of the reference genome listed in the GTF file; used to name STAR output files; ex.”M32”, “103”. | String |
| genome_build | Assembly accession (NCBI) or version (GENCODE) of the reference genome listed in the GTF file; used to name output files; ex. “GRCm39”, “GCF_003339765.1”.  | String |
| organism | Organism of the reference genome; used to name the output files; can be set to “Macaque”, “Mouse”, “Human”, or any other organism matching the reference genome. | String | 
| annotations_gtf | GTF file containing gene annotations; used to build the STAR reference files. | File |
| genome_fa | Genome FASTA file used for building indices. | File |
| biotypes | TSV file containing gene biotypes attributes to include in the modified GTF file; the first column contains the biotype and the second column contains “Y” to include or “N” to exclude the biotype; [GENCODE biotypes](https://www.gencodegenes.org/pages/biotypes.html) are used for GENCODE references and RefSeq biotypes are used for NCBI references. | File |

## BuildIndices tasks and tools

Overall, the BuildIndices workflow:
1. Checks inputs, modifies reference files, and creates STAR index.
2. Calculates chromosome sizes.
3. Builds reference bundle for bwa.

The tasks and tools used in the BuildIndices workflow are detailed in the table below. 

To see specific tool parameters, select the [workflow WDL link](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/build_indices/BuildIndices.wdl); then find the task and view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `docker: `.

| Task name | Tool | Software | Description | 
| --- | --- | --- | --- | 
| BuildStarSingleNucleus | [modify_gtf.py](https://github.com/broadinstitute/warp-tools/blob/develop/3rd-party-tools/build-indices/modify_gtf.py), STAR | [warp-tools](https://github.com/broadinstitute/warp-tools/tree/develop), [STAR](https://github.com/alexdobin/STAR) | Checks that the input GTF file contains input genome source, genome build version, and annotation version with correct build source information, modifies files for the STAR aligner, and creates STAR index file. |
| CalculateChromosomeSizes | faidx | [Samtools](http://www.htslib.org/) | Reads the genome FASTA file to create a FASTA index file that contains the genome chromosome sizes. |
| BuildBWAreference | index | [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) | Builds the reference bundle for the bwa aligner. |

#### 1. Check inputs, modify reference files, and create STAR index file

**Check inputs**

The BuildStarSingleNucleus task reads the input GTF file and verifies that the `genome_source`, `genome_build`, and `gtf_annotation_version` listed in the file match the input values provided to the pipeline.

**Modify reference files and create STAR index**

The BuildStarSingleNucleus task uses a custom python script, `[modify_gtf.py](https://github.com/broadinstitute/warp-tools/blob/develop/3rd-party-tools/build-indices/modify_gtf.py)`, and a list of biotypes ([example](https://github.com/broadinstitute/warp-tools/blob/develop/3rd-party-tools/build-indices/Biotypes.tsv)) to filter the input GTF file for only the biotypes indicated in the list with the value “Y” in the second column. The defaults in the custom code produce reference outputs that are similar to those built with 10x Genomics reference scripts.

The task uses the filtered GTF file and STAR `--runMode genomeGenerate` to generate the index file for the STAR aligner. Outputs of the task include the modified GTF and compressed STAR index files.

#### 2. Calculates chromosome sizes

The CalculateChromosomeSizes task uses Samtools to create and output a FASTA index file that contains the genome chromosome sizes, which can be used in downstream tools like SnapATAC2. 

#### 3. Builds reference bundle for bwa-mem2

The BuildBWAreference task uses the chromosome sizes file and bwa-mem2 to prepare the genome FASTA file for alignment and builds, compresses, and outputs the reference bundle for the bwa-mem2 aligner.

## Outputs

The following table lists the output variables and files produced by the pipeline.

| Output name | Filename, if applicable | Output format and description |
| ------ | ------ | ------ |
| snSS2_star_index | `modified_star2.7.10a-<organism>-<genome_source>-build-<genome_build>-<gtf_annotation_version>.tar` | TAR file containing a species-specific reference genome and GTF file for [STAR](https://github.com/alexdobin/STAR) alignment. |
| pipeline_version_out | `BuildIndices_v<pipeline_version>` | String describing the version of the BuildIndices pipeline used. |
| snSS2_annotation_gtf_modified | `modified_v<gtf_annotation_version>.annotation.gtf` | GTF file containing gene annotations filtered for selected biotypes. |
| reference_bundle | `bwa-mem2-2.2.1-<organism>-<genome_source>-build-<genome_build>.tar` | TAR file containing the reference index files for [BWA-mem](https://github.com/lh3/bwa) alignment. |
| chromosome_sizes | `chrom.sizes` | Text file containing chromosome sizes for the genome build. |

## Versioning and testing

All BuildIndices pipeline releases are documented in the [BuildIndices changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/build_indices/BuildIndices.changelog.md) and tested manually using [reference JSON files](https://github.com/broadinstitute/warp/tree/master/pipelines/skylab/build_indices).

## Consortia support
This pipeline is supported by the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and BICAN). 

If your organization also uses this pipeline, we would like to list you! Please reach out to us by contacting the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org).

## Feedback

Please help us make our tools better by contacting the [WARP Pipelines Team](mailto:warp-pipelines-help@broadinstitute.org) for pipeline-related suggestions or questions.