---
sidebar_position: 1
---

# Optimus Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [optimus_v5.1.1](https://github.com/broadinstitute/warp/releases?q=optimus&expanded=true) | November, 2021 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in warp or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

![Optimus_diagram](Optimus_diagram.png)

## Introduction to the Optimus workflow

Optimus is an open-source, cloud-optimized pipeline developed by the Data Coordination Platform (DCP) of the [Human Cell Atlas (HCA) Project](https://data.humancellatlas.org/) as well as the [BRAIN Initiative Cell Census Network](https://biccn.org/) (BICCN). It supports the processing of any 3' single-cell and single-nucleus count data generated with the [10x Genomics v2 or v3 assay](https://www.10xgenomics.com/solutions/single-cell/).

It is an alignment and transcriptome quantification pipeline that corrects cell barcodes (CBs), aligns reads to the genome, corrects Unique Molecular Identifiers (UMIs), generates a count matrix in a UMI-aware manner, calculates summary metrics for genes and cells, detects empty droplets, returns read outputs in BAM format, and returns cell gene counts in numpy matrix and Loom file formats.

In addition to providing commonly used metrics such as empty drop detection and mitochondrial reads, Optimus takes special care to **keep all reads in the output BAM that may be useful to the downstream user**, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development. 

Optimus has been validated for analyzing both human and mouse single-cell or single-nucleus datasets. It is currently optimized for samples less than 100 GB. Learn more in the [validation section](#validation-against-cell-ranger).

:::tip Want to use the Optimus pipeline for your publication?
Check out the [Optimus Publication Methods](./optimus.methods.md) to get started!
:::

## Quickstart table
The following table provides a quick glance at the Optimus pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | 10x Single Cell or Single Nucleus Expression (v2 and v3) | [10x Genomics](https://www.10xgenomics.com)
| Overall workflow  | Quality control module and transcriptome quantification module | Code available from [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence| GRCh38 human genome primary sequence and M21 (GRCm38.p6) mouse genome primary sequence | GENCODE [Human](https://www.gencodegenes.org/human/release_27.html) and [Mouse](https://www.gencodegenes.org/mouse/release_M21.html)
| Transcriptomic reference annotation | V27 GENCODE human transcriptome and M21 mouse transcriptome | GENCODE [Human](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz) and [Mouse](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz) |
| Aligner and transcript quantification | STARsolo (v.2.7.9a) | [Dobin, et al.,2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1) |
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which Optimus output is provided | [BAM](http://samtools.github.io/hts-specs/), Python numpy arrays (internal), Loom (generated with [Loompy v.3.0.6)](http://loompy.org/) |

## Set-up

### Optimus installation

To download the latest Opitmus release, see the release tags prefixed with "Optimus" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All Optimus pipeline releases are documented in the [Optimus changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/optimus/Optimus.changelog.md). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/develop/wreleaser).

If you’re running an Optimus workflow version prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the source code folder “website/pipelines/skylab/optimus).

Optimus can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. The Terra [Optimus Featured Workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline) contains the Optimus workflow, workflow configurations, required reference data and other inputs, and example testing data.


### Inputs

Optimus pipeline inputs are detailed in JSON format configuration files. There are five downsampled example configuration files available for running the pipeline:
*  [human_v2_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/human_v2_example.json): An example human 10x v2 single-cell dataset
*  [human_v3_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/human_v3_example.json): An example human 10x v3 single-cell dataset
*  [mouse_v2_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_example.json): An example mouse 10x v2 single-cell dataset
*  [mouse_v2_snRNA_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_snRNA_example.json): An example mouse v2 single-nucleus dataset

Additionally, there are multiple full-size example datasets available in the [test_optimus_full_datasets](https://github.com/broadinstitute/warp/tree/master/pipelines/skylab/optimus/example_inputs/test_optimus_full_datasets) folder. Unlike the example configuration files above, the full-size configuration files may not reflect updated Optimus parameters. However, you can still access the FASTQ files for each dataset at the Google bucket locations listed in the dataset configuration files.

#### Sample data input

Each 10x v2 and v3 3’ sequencing experiment generates triplets of FASTQ files for any given sample. Optimus takes the FASTQs list below as input. The pipeline is optimized for samples under 100 GB. To run larger samples, increasing the memory (the `machine_mem_mb` attribute) on the STARsoloFastq task.

1. Forward reads (`r1_fastq`) containing the unique molecular identifier (UMI) and cell barcode sequences
2. Reverse reads (`r2_fastq`) containing the alignable genomic information from the mRNA transcript
3. Optional index FASTQ (`i1_fastq`) containing the sample barcodes, when provided by the sequencing facility


:::tip Optimus is currently a single sample pipeline

However, it can take in multiple sets of FASTQs for a sample that has been split over multiple lanes of sequencing. For an example configuration file with multiple lanes, please see the [mouse_v2_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_example.json). Additionally, Optimus does not support demultiplexing even though it accepts index FASTQ files.
:::


#### Additional reference inputs

The example configuration files also contain metadata for the reference files, described in the  table below.

| Parameter name | Description | Optional strings (when applicable) |
| --- | --- | --- |
| whitelist | Cloud path to list of known CBs from [10x Genomics](https://www.10xgenomics.com/) that corresponds to the v2 or v3 chemistry. | NA |
| tar_star_reference | Cloud path to TAR file containing a species-specific reference genome and GTF; it is generated using the [BuildIndices workflow](https://github.com/broadinstitute/warp/tree/develop/pipelines/skylab/build_indices/BuildIndices.wdl). | NA |
| input_id | Unique identifier describing the biological sample or replicate that corresponds with the FASTQ files; can be a human-readable name or UUID. | NA |
| input_name | Optional string that can be used to further identify the original biological sample. | NA |
| input_id_metadata_field | Optional string describing, when applicable, the metadata field containing the input_id. | NA |
| input_name_metadata_field | Optional string describing, when applicable, the metadata field containing the input_name. | NA |
| annotations_gtf | Cloud path to GTF containing gene annotations used for gene tagging (must match GTF in STAR reference). | NA |
| chemistry | Optional string description of whether data was generated with 10x v2 or v3 chemistry. Optimus validates this string. If the string does not match one of the optional strings, the pipeline will fail. You can remove the checks by setting "force_no_check = true" in the input JSON | "tenX_v2" (default) or "tenX_v3". |
| counting_mode | String description of whether data is single-cell or single-nucleus | "sc_rna" or "sn_rna". |
| output_bam_basename | Optional string used for the output BAM file basename; the default is input_id. | NA |
| use_strand_info | Optional string for reading stranded data. Default is "false"; set to "true" to count reads in stranded mode. | "true" or "false" (default) |
| emptydrops_lower | UMI threshold for emptyDrops detection; default is 100. | NA |

#### Pseudogene handling
Optimus reference files are downloaded directly from GENCODE (see Quickstart table) and are not modified to remove pseudogenes. This is in contrast to the [references created for Cell Ranger](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/release-notes/references#header) which remove pseudogenes and small RNAs.

In the case of multi-mapped pseudogenes, Optimus and Cell Ranger will produce different results. Optimus does not count multi-mapped reads in the final count matrix, whereas Cell Ranger will keep potential multi-mapped reads because it does not identify the pseudogene reads.

#### Sample inputs for analyses in a Terra Workspace

The Optimus pipeline is currently available on the cloud-based platform Terra. After registering, you can access the Featured Workspace using this address: [https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline). The workspace is preloaded with instructions and sample data. Please view the [Support Center](https://support.terra.bio/hc/en-us) for more information on using the Terra platform.

## Optimus tasks and tools

The [Optimus workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) imports individual "tasks," also written in  WDL script, from the WARP [tasks folder](https://github.com/broadinstitute/warp/blob/master/tasks/skylab). 

Overall, the Optimus workflow:
1. Corrects CBs, aligns reads, corrects UMIs, and counts genes with STAR v.2.7.9a.
1. Calculates gene metrics.
1. Calculates cell metrics.
1. Converts the Star Output into NPY and NPZ arrays.
1. Runs emptyDrops.
1. Merges gene counts, metrics, and emptyDrops data into a Loom-formatted matrix.


The tools each Optimus task employs are detailed in the table below. 

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `String docker =`.

| Task name and WDL link | Tool | Software | Description | 
| --- | --- | --- | ------------------------------------ | 
| [StarAlign.STARsoloFastq (alias = STARsoloFastq)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlign.wdl) | STAR | [Star v2.7.9a](https://github.com/alexdobin/STAR) | Uses the input FASTQ files to perform cell barcode correction, adaptor trimming, alignment, gene annotation, UMI correction, and gene counting. |
| [Metrics.CalculateGeneMetrics (alias = GeneMetrics)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) | TagSort | sctools | Sorts the BAM file by gene using the cell barcode (CB), molecule barcode (UB) and gene ID (GX) tags and computes gene metrics. | 
| [Metrics.CalculateCellMetrics (alias = CellMetrics)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) | TagSort | sctools | Sorts the BAM file by cell using the cell barcode (CB), molecule barcode (UB) and gene ID (GX) tags and computes cell metrics. |
| [StarAlign.ConvertStarOutput (alias = ConvertOutputs)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlign.wdl) | create-npz-output.py | Python3 | Creates a compressed raw NPY or NPZ file containing the STARsolo output features (NPY), barcodes (NPZ) and counts (NPZ). | 
| [RunEmptyDrops.RunEmptyDrops](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/RunEmptyDrops.wdl) | npz2rds.sh, emptyDropsWrapper.R, emptyDrops | [DroploetUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) | Runs custom scripts to convert the NPY and NPZ files to RDS and then uses emptyDrops to identify empty lipid droplets. |
|  [LoomUtils.OptimusLoomGeneration](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/LoomUtils.wdl) | create_loom_optimus.py | Python3 | Merges the gene counts, cell metrics, gene metrics, and emptyDrops data into a Loom formatted cell-by-gene matrix. |


More information about the different tags used to flag the data can be found in the [Bam_tags documentation](./Bam_tags.md).


#### 1. CB correction, read trimming, alignment, gene annotation, UMI correction, and gene counting

The [STARsoloFastq task](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlign.wdl) task uses the array of forward and reverse read FASTQ files to perform cell barcode correction, alignment, gene annotation, and counting. 

**CB correction**

Although the function of the CBs is to identify unique cells, barcode errors can arise during sequencing, such as the incorporation of the barcode into contaminating DNA or sequencing and PCR errors. This makes it difficult to distinguish unique cells from artifactual appearances of the barcode. The STARsoloFastq task uses the STAR aligner to evaluate barcode errors by comparing the R1 FASTQ sequences against a 10x chemistry-specific whitelist of known barcode sequences. 

Corrected barcodes are those that come within one edit distance ([Hamming  distance](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410656/)) of matching the whitelist of barcode sequences. This is specified in the STAR parameter  `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`, which allows multiple matches in the whitelist with 1 mismatched base and uses posterior probability calculation to choose one of the matches. Additionally, it adds pseudocounts of 1 to all whitelist barcodes and allows multi-matching of CBs with N-bases to the whitelist.

Correct barcodes are assigned a “CB” tag in the BAM. Uncorrectable barcodes (with more than one error) are preserved and given a “CR” (Cell barcode Raw) tag. Cell barcode quality scores are also preserved in the file under the “CY” tag. 

**Read trimming**

Read trimming removes Illumina adaptor sequences. This is set to match the read trimming performed by CellRanger4 and is specified using the parameter `--clipAdapterType CellRanger4` and  `--outFilterScoreMin 30`.

**Alignment**

STAR maps barcoded reads to the genome primary assembly reference (see the Quickstart table above for version information). This reference is generated using the [BuildIndices workflow](https://github.com/broadinstitute/warp/tree/develop/pipelines/skylab/build_indices/BuildIndices.wdl). The strandedness for alignment is specified with the `--soloStrand` parameter, which is set to unstranded by default. 

**Gene annotation**

Prior to gene counting, STARsolo adds gene annotations which will vary depending on the counting_mode ("sc_rna" or "sn_rna") specified in the workflow. With sc_rna, STARsolo runs with the “Gene” COUNTING_MODE, whereas with sn_rna, it runs the “GeneFull” COUNTING_MODE.

Genes that overlap an alignment are stored with the GX BAM tag; for sc_rna mode, this will include the gene that corresponds to an exon or UTR, whereas for sn_rna mode, this will include the gene corresponding to an exon, UTR, and intron.

All tags are detailed in the pipeline's [BAM_tag documentation](./Bam_tags.md).

**UMI correction and gene counting**

UMIs are designed to distinguish unique transcripts present in the cell at lysis from those arising from PCR amplification of these same transcripts. But, like CBs, UMIs can also be incorrectly sequenced or amplified. 

By specifying the --soloUMIdedup 1MM_Directional_UMItools, STARsolo applies a network-based, "directional" correction method ([Smith, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/)) to account for such errors. Additionally, the UMI correction task requires the length of the UMI, which differs between v2 and v3 chemistry; v2 is 10 bp whereas v3 is 12 bp. The task will add a 'UB' tag for UMI-corrected barcodes. 

Deduplicated UMIs are counted towards their assigned gene/cells, producing a raw count matrix.


**STARsolo outputs**

The task’s output includes a coordinate-sorted BAM file containing the cell barcode-corrected reads and SAM attributes UB UR UY CR CB CY NH GX GN. Additionally, after counting, the task outputs three intermediate TSV files (features, barcodes, and matrix) used for downstream empty droplet detection and Loom matrix generation.   

#### 2. Calculate gene metrics

The [CalculateGeneMetrics](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) task uses [sctools](https://github.com/HumanCellAtlas/sctools) to calculate summary metrics that help assess the quality of the data output each time this pipeline is run. These metrics are included in the output Loom matrix. A detailed list of these metrics is found in the [Optimus Count Matrix Overview](./Loom_schema.md).

#### 3. Calculate cell metrics

The [CalculateCellMetrics](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) task uses [sctools](https://github.com/HumanCellAtlas/sctools) to calculate summary metrics that help assess the per-cell quality of the data output each time this pipeline is run. These metrics are included in the output Loom matrix. A detailed list of these metrics is found in the [Optimus Count Matrix Overview](./Loom_schema.md).


#### 4. Convert STAR output

The [ConvertStarOutput](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlign.wdl) task uses a custom python script to convert the matrix, features, and barcodes output from STARsolo into an NPY (features and barcodes)- and NPZ (the matrix)-formatted file for downstream empty drops detection and Loom matrix generation. 

#### 5. Run emptyDrops

Empty droplets are lipid droplets that did not encapsulate a cell during 10x sequencing, but instead acquired cell-free RNA (secreted RNA or RNA released during cell lysis) from the solution in which the cells resided ([Lun, et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/?term=30902100). This ambient RNA can serve as a substrate for reverse transcription, leading to a small number of background reads. The Optimus pipeline calls the [RunEmptyDrops](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/RunEmptyDrops.wdl) task which uses the [DropletUtils v1.2.1](http://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package to flag CBs that represent empty droplets rather than cells. A cell will be flagged if it contains fewer than 100 molecules.
emptyDrops metrics for single-cell data (not single-nucleus; see note below) are stored in the columns of the output Loom matrix.  Read full details for all the metrics in the [Optimus Matrix Overview](./Loom_schema.md).

EmptyDrops relies on a visual knee point inflection (described in [Lun et al. (2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)) to differentiate ambient-like cells from empty droplets. **If the single-cell data (counting mode set to `sc_rna`) does not produce a knee point inflection when running emptyDrops, the Loom columns for emptyDrops data will contain "NA"s.**

:::warning RunEmptyDrops output not included for single-nucleus data
Often snRNAseq data does not produce a visual knee point inflection when running emptyDrops and the tool can not accurately distinguish ambient-like cells from empty droplets. For this reason, emptyDrops is not used if Optimus counting_mode is set to `sn_rna`, and the output Loom matrix will not contain emptyDrops metrics.
:::


#### 5.  Matrix construction

The [OptimusLoomGeneration](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/LoomUtils.wdl) task uses a custom python script to merge the converted STARsolo count matrix, the emptyDrops results, and the cell and gene metrics into a Loom-formatted cell-by-gene matrix. Read full details for all the metrics in the [Optimus Count Matrix Overview](./Loom_schema.md).



#### 9. Outputs

Output files of the pipeline include:

1. Cell x Gene unnormalized, but UMI-corrected, count matrices in Loom format.
2. Unfiltered, sorted BAM file with barcode and downstream analysis [tags](./Bam_tags.md).
3. Cell metadata, including cell metrics.
4. Gene metadata, including gene metrics.

The following table lists the output files produced from the pipeline. For samples that have sequenced over multiple lanes, the pipeline will output one merged version of each listed file.

| Output Name | Filename, if applicable | Output Type |Output Format |
| ------ |------ | ------ | ------ |
| pipeline_version | | Version of the processing pipeline run on this data. | String |
| bam | `<input_id>.bam` | Aligned BAM | BAM |
| matrix | Converted sparse matrix file from the ConvertStarOutput task. | NPZ |
| matrix_row_index | sparse_counts_row_index.npy | Index of cells in count matrix. | NPY |
| matrix_col_index | sparse_counts_col_index.npy | Index of genes in count matrix. | NPY |
| cell_metrics | cell-metrics.csv.gz | cell metrics | compressed csv | Matrix of metrics by cells. |
| gene_metrics | gene-metrics.csv.gz | gene metrics | compressed csv | Matrix of metrics by genes. |
| cell_calls | empty_drops_result.csv | emptyDrops results from the RunEmptyDrops task. | CSV |
| loom_output_file | `<input_id>.loom` | Loom | Loom | Loom file with count data and metadata. | N/A |

The Loom matrix is the default output. See the [create_loom_optimus.py](https://github.com/broadinstitute/warp/blob/master/dockers/skylab/loom-output/create_loom_optimus.py) for the detailed code. This matrix contains the unnormalized (unfiltered), UMI-corrected count matrices, as well as the gene and cell metrics detailed in the [Optimus Count Matrix Overview](./Loom_schema.md).

The matrix is compatible with multiple downstream community analysis tools. For a tutorial on using the Optimus matrix with [Seurat](https://satijalab.org/seurat/index.html), [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Cumulus](https://cumulus.readthedocs.io/en/latest/index.html), or [Pegasus](https://pegasus.readthedocs.io/en/stable/#), see the public [Intro-to-HCA-data-on-Terra workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/Intro-to-HCA-data-on-Terra) (login required) and its accompanying [step-by-step guide](https://support.terra.bio/hc/en-us/articles/360060041772).


## Validation against Cell Ranger
Optimus has been validated for processing both human and mouse single-cell and single-nucleus data (see links to validation reports in the table below). For each validation, Optimus results are compared to those of Cell Ranger (see the [FAQ](#faqs) for more on Cell Ranger comparisons).

| Workflow configuration | Link to Report |
| --- | --- |
| Human 10x v2 single-cell | [Report](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/benchmarking/v1_Apr2019/optimus_report.rst)|
| Mouse 10x v2 single-cell | [Report](https://docs.google.com/document/d/1_3oO0ZQSrwEoe6D3GgKdSmAQ9qkzH_7wrE7x6_deL10/edit) |
| Human and mouse 10x v3 single-cell | [Report](https://docs.google.com/document/d/1-hwfXkqtL8MblgDWFzk-HsVRYiy4PS8ZhJqAGlHBWYE/edit#heading=h.4uokn64v1s5m)
|Human and Mouse 10x v2/v3 single-nucleus | [Report](https://docs.google.com/document/d/1rv2M7vfpOzIOsMnMfNyKB4HV18lQ9dnOGHK2tPikiH0/edit) |
| Optimus STARsolo (v5.0.0 and later) | [Report](https://docs.google.com/document/d/1B6Ux6HICD4ZL4Z0TG9LO-X43gOdF3sbq5qw_L2GA6fg/edit) |


## Versioning

All Optimus pipeline releases are documented in the [Optimus changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/optimus/Optimus.changelog.md).


## Citing the Optimus Pipeline
Please identify the pipeline in your methods section using the Optimus Pipeline's [SciCrunch resource identifier](https://scicrunch.org/scicrunch/Resources/record/nlx_144509-1/SCR_018908/resolver?q=SCR_018908&l=SCR_018908).
* Ex: *Optimus Pipeline (RRID:SCR_018908)*

## Consortia Support
This pipeline is supported and used by the [Human Cell Atlas](https://www.humancellatlas.org/) (HCA) project and the[BRAIN Initiative Cell Census Network](https://biccn.org/) (BICCN). 

If your organization also uses this pipeline, we would like to list you! Please reach out to us by contacting [Kylee Degatano](mailto:kdegatano@broadinstitute.org).

## Feedback

Please help us make our tools better by contacting [Kylee Degatano](mailto:kdegatano@broadinstitute.org) for pipeline-related suggestions or questions.


## FAQs

:::note Question Can I run Optimus in Terra?

Yes! We have a Terra workspace that is preconfigured with the latest Optimus workflow and is preloaded with human and mouse sample data. You can access the [workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline). You will need a Google account to set up Terra. Please see [Terra Support](https://support.terra.bio/hc/en-us) for documents on getting started.
:::

:::note Question Is the output count matrix filtered or normalized?

No, we do not filter. We keep as much data as possible so that the researcher can make their own filtering and normalization choices. We do, however, output some information that may be helpful for filtering, like UMI counts per cell and calls on whether or not a cell is empty from emptyDrops software. For the emptyDrops call, a cell will be flagged as possibly empty if it contains fewer than 100 molecules.
:::

:::note Question How does the workflow change when using the single-cell RNA-seq (counting_mode = 'sc_rna') vs. the single-nucleus (counting_mode = 'sn_rna') parameters?

The counting_mode parameter is used to specify the STARsolo COUNTING_MODE; when sn_rna is specified, STARsolo will tag gene exons, UTRs, AND introns with the GX tag.  Additionally, the Optimus uses the counting_mode to determine whether to run emptyDrops; no emptyDrops data is calculated for the sn_rna mode.
:::

:::note Question Where can I find example Optimus datasets and parameters to test the pipeline?

There are four example configuration JSON files available for you to test the pipeline- the [human_v2_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/human_v2_example.json), the [human_v3_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/human_v3_example.json), the [mouse_v2_snRNA_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_snRNA_example.json), and the [mouse_v2_snRNA_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_snRNA_example.json)(see the Inputs section). Each of these configuration files can be run in the Optimus Featured Workspace in Terra at https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline, but you should note that the workspace comes preloaded with the same data and configurations. We also have multiple example datasets available in the [test_optimus_full_datasets](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/test_optimus_full_datasets/) folder. These datasets contain configuration files listing the cloud location of the dataset FASTQ files; however, the configuration files may not be updated with all the workflow parameters. For the most up-to-date configuration examples, see the four example files listed above.
:::

:::note Question What outputs are expected if my sample has been sequenced over multiple lanes?

The Optimus pipeline is a single sample pipeline, but it can accept multiple FASTQ files if a sample is sequenced across lanes. In this case, the pipeline will merge the results from each lane into single output files. There will only be one merged file for each output type (i.e one Loom, etc.). If you would like to view an example configuration file for a multi-lane dataset, please see the [mouse_v2_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_example.json).  Additionally, you can view sample outputs in the Optimus featured workspace on [Terra](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline).
:::

:::note Question How do I find which parameters and Docker images were used for the different tasks (i.e. STAR alignment, emptyDrops, etc.)

Parameters are listed in each task WDL. For a list of the tasks, see the table in the [Task Summary Section](#optimus-task-summary). Select the link for the task of interest and then view the parameters in the task WDL "command {}" section. For the task Docker image, see task WDL "# runtime values" section; the Docker is listed as "String docker =  ". If you want to learn more about all the different parameters available for a software tool, please select the relevant link in the table's "Tool" column.
:::

:::note Question Does Optimus have any read length requirements?
For Read 1 sequences, the only minimum requirement is that reads are the combined lengths of the cell barcode and UMIs (which will vary between 10x V1, V2, and V3 chemistry).

For Read 2 sequences, there is no read length requirement and read lengths will vary.
:::

:::note Question How does Optimus compare to Cell Ranger?

Cell Ranger is a commonly used set of analysis pipelines developed by [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). Optimus and Cell Ranger share many features and additionally, Optimus results are validated against Cell Ranger results (see our [human validation report](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/benchmarking/v1_Apr2019/optimus_report.rst)). 

*So why develop an independent pipeline for 10x data analyses?*

For three reasons:
1) Need for an open-source, cloud-optimized pipeline. When Optimus was developed, Cell Ranger software was not yet open source, nor was it optimized for the cloud. To date, the Cell Ranger open-source code is still not regularly updated with Cell Ranger releases. In consequence, using the latest Cell Ranger (which is not open source yet) limits our ability to harness the breadth of tools available in the scientific community.

2) Flexibility to process data similar, but not identical, to 10x. We wanted the ability to evolve our pipeline to process non-10x data types that might use similar features such as combinatorial indexing.

3) Addition of metrics. We wanted the pipeline to calculate key metrics that would be useful to the scientific community, such as emptyDrops calculations, mitochondrial read metrics, etc.

*Reference differences between Optimus and Cell Ranger*

Unlike Cell Ranger references, Optimus references are downloaded directly from GENCODE and not modified to remove pseudogenes and small RNAs. Learn more about Cell Ranger references on the [10x website](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/release-notes/references#header). 

In the case of multi-mapped pseudogenes, Optimus and Cell Ranger will produce different results. Optimus does not count multi-mapped reads in the final count matrix, whereas Cell Ranger will keep potential multi-mapped reads because it does not identify the pseudogene reads.
:::



