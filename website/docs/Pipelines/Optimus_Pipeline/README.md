---
sidebar_position: 1
---

# Optimus Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [optimus_v5.7.0](https://github.com/broadinstitute/warp/releases?q=optimus&expanded=true) | February, 2023 | Elizabeth Kiernan | Please file GitHub issues in warp or contact [the WARP team](mailto:warp-pipelines-help@broadinstitute.org) |

![Optimus_diagram](Optimus_diagram.png)

## Introduction to the Optimus workflow

Optimus is an open-source, cloud-optimized pipeline developed by the Data Coordination Platform (DCP) of the [Human Cell Atlas (HCA) Project](https://data.humancellatlas.org/) as well as the [BRAIN Initiative Cell Census Network](https://biccn.org/) (BICCN). It supports the processing of any 3' single-cell and single-nucleus count data generated with the [10x Genomics v2 or v3 assay](https://www.10xgenomics.com/solutions/single-cell/).

It is an alignment and transcriptome quantification pipeline that corrects cell barcodes (CBs), aligns reads to the genome, corrects Unique Molecular Identifiers (UMIs), generates a count matrix in a UMI-aware manner, calculates summary metrics for genes and cells, detects empty droplets, returns read outputs in BAM format, and returns cell gene counts in numpy matrix and Loom file formats.

In addition to providing commonly used metrics such as empty drop detection and mitochondrial reads, Optimus takes special care to **keep all reads in the output BAM that may be useful to the downstream user**, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development. 

Optimus has been validated for analyzing both human and mouse single-cell or single-nucleus datasets. Learn more in the [validation section](#validation-against-cell-ranger).

:::tip Want to use the Optimus pipeline for your publication?
Check out the [Optimus Publication Methods](./optimus.methods.md) to get started!
:::

## Quickstart table
The following table provides a quick glance at the Optimus pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | 10x single cell or single nucleus expression (v2 and v3) | [10x Genomics](https://www.10xgenomics.com)
| Overall workflow  | Quality control module and transcriptome quantification module | Code available from [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence and M21 (GRCm38.p6) mouse genome primary sequence | GENCODE [human reference files](https://www.gencodegenes.org/human/release_27.html) and [mouse reference files](https://www.gencodegenes.org/mouse/release_M21.html)
| Transcriptomic reference annotation | V27 GENCODE human transcriptome and M21 mouse transcriptome | GENCODE [human GTF](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz) and [mouse GTF](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz) |
| Aligner and transcript quantification | STARsolo | [Dobin, et al.,2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1) |
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which Optimus output is provided | [BAM](http://samtools.github.io/hts-specs/), Python numpy arrays (internal), Loom (generated with [Loompy v.3.0.6)](http://loompy.org/) |

## Set-up

### Optimus installation

To download the latest Optimus release, see the release tags prefixed with "Optimus" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All Optimus pipeline releases are documented in the [Optimus changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/optimus/Optimus.changelog.md). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

If you’re running an Optimus workflow version prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the source code folder “website/pipelines/skylab/optimus").

Optimus can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. The Terra [Optimus Featured Workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline) contains the Optimus workflow, workflow configurations, required reference data and other inputs, and example testing data.


### Inputs

Optimus pipeline inputs are detailed in JSON format configuration files. There are five downsampled example configuration files available for running the pipeline:
*  [human_v2_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/human_v2_example.json): An example human 10x v2 single-cell dataset
*  [human_v3_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/human_v3_example.json): An example human 10x v3 single-cell dataset
*  [mouse_v2_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_example.json): An example mouse 10x v2 single-cell dataset
*  [mouse_v2_snRNA_example](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_snRNA_example.json): An example mouse v2 single-nucleus dataset

Additionally, there are multiple full-size example datasets available in the [test_optimus_full_datasets](https://github.com/broadinstitute/warp/tree/master/pipelines/skylab/optimus/example_inputs/test_optimus_full_datasets) folder. Unlike the example configuration files above, the full-size configuration files may not reflect updated Optimus parameters. However, you can still access the FASTQ files for each dataset at the Google bucket locations listed in the dataset configuration files.

#### Sample data input

Each 10x v2 and v3 3’ sequencing experiment generates triplets of FASTQ files for any given sample. Optimus takes the FASTQs listed below as input.

1. Forward reads (`r1_fastq`) containing the unique molecular identifier (UMI) and CB sequences
2. Reverse reads (`r2_fastq`) containing the alignable genomic information from the mRNA transcript
3. Optional index FASTQ (`i1_fastq`) containing the sample barcodes, when provided by the sequencing facility


:::tip Optimus is currently a single sample pipeline

However, it can take in multiple sets of FASTQs for a sample that has been split over multiple lanes of sequencing. For an example configuration file with multiple lanes, please see the [mouse_v2_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/example_inputs/mouse_v2_example.json). Additionally, Optimus does not support demultiplexing even though it accepts index FASTQ files.
:::


#### Additional reference inputs

The example configuration files also contain metadata for the reference files, described in the table below.

| Parameter name | Description | Optional attributes (when applicable) |
| --- | --- | --- |
| whitelist |  List of known CBs; the workflow automatically selects the [10x Genomics](https://www.10xgenomics.com/) whitelist that corresponds to the v2 or v3 chemistry based on the input `tenx_chemistry_version`. A custom whitelist can also be provided if the input data was generated with a chemistry different from 10x Genomics v2 or v3. To use a custom whitelist, set the input `ignore_r1_read_length` to "true". | N/A |
| tar_star_reference | TAR file containing a species-specific reference genome and GTF; it is generated using the [BuildIndices workflow](https://github.com/broadinstitute/warp/tree/develop/pipelines/skylab/build_indices/BuildIndices.wdl). | N/A |
| input_id | Unique identifier describing the biological sample or replicate that corresponds with the FASTQ files; can be a human-readable name or UUID. | N/A |
| input_name | Optional string that can be used to further identify the original biological sample. | N/A |
| input_id_metadata_field | Optional string describing, when applicable, the metadata field containing the input_id. | N/A |
| input_name_metadata_field | Optional string describing, when applicable, the metadata field containing the input_name. | N/A |
| annotations_gtf | GTF containing gene annotations used for gene tagging (must match GTF in STAR reference). | N/A |
| tenx_chemistry_version | Integer that specifies if data was generated with 10x v2 or v3 chemistry. Optimus validates this chemistry by examining the UMIs and CBs in the first read 1 FASTQ file. If the chemistry does not match, the pipeline will fail. You can remove the check by setting "ignore_r1_read_length = true" in the input JSON. | 2 or 3 |
| mt_genes | Optional file containing mitochondrial gene names for a specific species. This is used for calculating gene metrics. | N/A |
| counting_mode | String describing whether data is single-cell or single-nucleus. Single-cell mode counts reads aligned to the gene transcript, whereas single-nucleus counts whole transcript to account for nuclear pre-mRNA. | "sc_rna" or "sn_rna" |
| output_bam_basename | String used as a basename for output BAM file; the default is set to the string used for the `input_id` parameter. | N/A |
| use_strand_info | Optional string for reading stranded data. Default is "false"; set to "true" to count reads in stranded mode. | "true" or "false" (default) |
| ignore_r1_read_length | Boolean that overrides a check on the 10x chemistry. Default is set to false. If true, the workflow will not ensure that the 10x_chemistry_version input matches the chemistry in the read 1 FASTQ. | "true" or "false" (default) | 
| emptydrops_lower | UMI threshold for emptyDrops detection; default is 100. | N/A |
| count_exons | Boolean indicating if the workflow should calculate exon counts **when in single-nucleus (sn_rna) mode**. If true, this option will output an additional layer for the Loom file. By default, it is set to "false". If the parameter is true and used with sc_rnamode, the workflow will return an error. | "true" or "false" (default) |

#### Pseudogene handling
The example Optimus reference files are downloaded directly from GENCODE (see Quickstart table) and are not modified to remove pseudogenes. This is in contrast to the [references created for Cell Ranger](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/release-notes/references#header) which remove pseudogenes and small RNAs.

In the case of multi-mapped pseudogenes, Optimus and Cell Ranger will produce different results. Optimus does not count multi-mapped reads in the final count matrix, whereas Cell Ranger will keep potential multi-mapped reads because it does not identify the pseudogene reads.

#### Sample inputs for analyses in a Terra Workspace

The Optimus pipeline is currently available on the cloud-based platform Terra. After registering, you can access the [Optimus featured workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline). The workspace is preloaded with instructions and sample data. Please view the [Support Center](https://support.terra.bio/hc/en-us) for more information on using the Terra platform.

## Optimus tasks and tools

The [Optimus workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) imports individual "tasks," also written in  WDL script, from the WARP [tasks folder](https://github.com/broadinstitute/warp/blob/master/tasks/skylab). 

Overall, the Optimus workflow:
1. Checks inputs 
1. Partitions FASTQs by CB.
1. Corrects CBs, aligns reads, corrects UMIs, and counts genes with STAR.
1. Merges the Star outputs into NPY and NPZ arrays.
1. Calculates gene metrics.
1. Calculates cell metrics.
1. Runs emptyDrops.
1. Merges gene counts, metrics, and emptyDrops data into a Loom-formatted matrix.


The tools each Optimus task employs are detailed in the table below. 

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `String docker =`.

| Task name and WDL link | Tool | Software | Description | 
| --- | --- | --- | ------------------------------------ | 
| [OptimusInputChecks.checkOptimusInput](https://github.com/broadinstitute/warp/blob/develop/tasks/skylab/CheckInputs.wdl) | Custom script | Bash | Validates the tenx_chemistry_version and counting_mode inputs. For the tenx_chemistry_version, the script verifies that the length of the CBs and UMIs in the first read1 FASTQ file matches the length expected by 10x Genomics v2 and v3 chemistry. It then selects and outputs the appropriate whitelist to use for analysis in the rest of the workflow. |
| [FastqProcessing.FastqProcessing (alias = SplitFastq)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/FastqProcessing.wdl) | fastqprocess | [warp-tools](https://github.com/broadinstitute/warp-tools) | Partitions the input FASTQ files by CB to create an array of FASTQ files that are each ~ 30 GB. This task is skipped if the input files are smaller than 30 GB. | 
| [StarAlign.STARsoloFastq (alias = STARsoloFastq)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlign.wdl) | STAR | [Star](https://github.com/alexdobin/STAR) | Uses the partitioned FASTQ files to perform CB correction, adaptor trimming, alignment, gene annotation, UMI correction, and gene counting. |
| [Merge.MergeSortBamFiles (alias= MergeBam)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/MergeSortBam.wdl) | MergeSamFiles | Picard | Merges the array of BAM files into a single BAM and sorts in coordinate order. |
| [StarAlign.MergeStarOutput (alias = MergeStarOutputs)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlign.wdl) | create-npz-output.py | Python3 | Creates a compressed raw NPY or NPZ file containing the STARsolo output features (NPY), barcodes (NPZ) and counts (NPZ). | 
| [Metrics.CalculateGeneMetrics (alias = GeneMetrics)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) | TagSort | [warp-tools](https://github.com/broadinstitute/warp-tools) | Sorts the BAM file by gene using the cell barcode (CB), molecule barcode (UB) and gene ID (GX) tags and computes gene metrics. | 
| [Metrics.CalculateCellMetrics (alias = CellMetrics)](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) | TagSort | [warp-tools](https://github.com/broadinstitute/warp-tools) | Sorts the BAM file by cell using the cell barcode (CB), molecule barcode (UB) and gene ID (GX) tags and computes cell metrics. |
| [RunEmptyDrops.RunEmptyDrops](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/RunEmptyDrops.wdl) | npz2rds.sh, emptyDropsWrapper.R, emptyDrops | [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) | Runs custom scripts to convert the NPY and NPZ files to RDS and then uses emptyDrops to identify empty lipid droplets. This step only runs when `counting_mode` = "sc_rna".|
|  [LoomUtils.OptimusLoomGeneration](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/LoomUtils.wdl) | create_loom_optimus.py | Python3 | Merges the gene counts, cell metrics, gene metrics, and emptyDrops data into a Loom formatted cell-by-gene matrix. The Loom contains exon counts when using sc_rna mode, and whole-gene counts when running in sn_rna mode. It optionally contains an additional layer for exon counts when running sn_rna mode with `exon_counts` set to true. |


More information about the different tags used to flag the data can be found in the [Bam_tags documentation](./Bam_tags.md).

#### 1. Checks inputs
The [checkOptimusInput](https://github.com/broadinstitute/warp/blob/develop/tasks/skylab/CheckInputs.wdl) task validates the `tenx_chemistry_version` and `counting_mode` inputs to ensure that the attributes are selected. For the `tenx_chemistry_version`, the script verifies that the length of the CBs and UMIs in the first read1 FASTQ file matches the length expected by 10x Genomics v2 and v3 chemistry. It then selects and outputs the appropriate whitelist to use for analysis in the rest of the workflow. For the `counting_mode` input, the task verifies that the input is set to either the string "sc_rna" or "sn_rna". If these checks fail, the workflow will fail. 

#### 2. Partition CBs
To optimally scale for large files, the workflow first uses the [FastqProcessing task](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/FastqProcessing.wdl) to partition the input array of forward and reverse read FASTQ files by uncorrected CBs. The resulting array contains FASTQ files that are \~30 GB in size. If the input files are smaller than 30 GB, this task is omitted.

#### 3. Correct CBs, trims reads, align, annotate genes, correct UMIs, and count genes
The [STARsoloFastq task](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlign.wdl) uses the partitioned FASTQ files to perform CB correction, alignment, gene annotation, and counting. 

**CB correction**

Although the function of the CBs is to identify unique cells, barcode errors can arise during sequencing, such as the incorporation of the barcode into contaminating DNA or sequencing and PCR errors. This makes it difficult to distinguish unique cells from barcode artifacts. The STARsoloFastq task uses the STAR aligner to evaluate barcode errors by comparing the R1 FASTQ sequences against a 10x chemistry-specific whitelist of known barcode sequences. 

Corrected barcodes are those that come within one edit distance ([Hamming  distance](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410656/)) of matching the whitelist of barcode sequences. This is specified in the STAR parameter  `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`, which allows multiple matches in the whitelist with 1 mismatched base and uses posterior probability calculation to choose one of the matches. Additionally, it adds pseudocounts of 1 to all whitelist barcodes and allows multi-matching of CBs with N-bases to the whitelist.

Correct barcodes are assigned a “CB” tag in the BAM. Uncorrectable barcodes (with more than one error) are preserved and given a “CR” (Cell barcode Raw) tag. CB quality scores are also preserved in the file under the “CY” tag. 

**Read trimming**

Read trimming removes Illumina adapter sequences. This is set to match the read trimming performed by CellRanger4 and is specified using the parameter `--clipAdapterType CellRanger4` and  `--outFilterScoreMin 30`.

**Alignment**

STAR maps barcoded reads to the genome primary assembly reference (see the Quickstart table above for version information). The example references for Optimus were generated using the [BuildIndices workflow](https://github.com/broadinstitute/warp/tree/develop/pipelines/skylab/build_indices/BuildIndices.wdl). The strandedness for alignment is specified in STAR with the `--soloStrand` parameter, which is set to unstranded by default. 

**Gene annotation**

Prior to gene counting, STARsolo adds gene annotations which will vary depending on the counting_mode ("sc_rna" or "sn_rna") specified in the Optimus workflow. With `sc_rna`, STARsolo runs with the “Gene” COUNTING_MODE, which is specific to exons.  With the `sn_rna` mode, STARsolo runs the “GeneFull” COUNTING_MODE and has the additional option to run the "Gene" mode when the input parameter `count_exons` is set to "true". 

Genes that overlap an alignment are stored with the GX BAM tag; for sc_rna mode, this will include the gene that corresponds to an exon or UTR, whereas for sn_rna mode, this will include the gene corresponding to an exon, UTR, and intron.

All tags are detailed in the pipeline's [BAM_tag documentation](./Bam_tags.md).

The resulting BAM files are merged together into a single BAM using the [MergeSortBamFiles (alias= MergeBam)](https://github.com/broadinstitute/warp/blob/develop/tasks/skylab/MergeSortBam.wdl).

**UMI correction and gene counting**

UMIs are designed to distinguish unique transcripts present in the cell at lysis from those arising from PCR amplification of these same transcripts. But, like CBs, UMIs can also be incorrectly sequenced or amplified. 

By specifying the `--soloUMIdedup 1MM_Directional_UMItools`, STARsolo applies a network-based, "directional" correction method ([Smith, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/)) to account for such errors. Additionally, the UMI correction task requires the length of the UMI, which differs between v2 and v3 chemistry; v2 is 10 bp whereas v3 is 12 bp. The task will add a 'UB' tag for UMI-corrected barcodes. 

Deduplicated UMIs are counted towards their assigned gene/cells, producing a raw count matrix. 

#### 4. Merge STARsolo count matrices

The STARsolo output will include a features, barcodes, and matrix TSV for each of the partitioned FASTQ input. The [MergeStarOutput task](https://github.com/broadinstitute/warp/blob/develop/tasks/skylab/StarAlign.wdl) merges each respective TSV. It uses a custom python script to convert the merged matrix, features, and barcodes output from STARsolo into an NPY (features and barcodes)- and NPZ (the matrix)-formatted file for downstream empty drops detection and Loom matrix generation. 

**STARsolo outputs**

The task’s output includes a merged, coordinate-sorted BAM file containing the CB-corrected reads and SAM attributes UB UR UY CR CB CY NH GX GN. Additionally, after counting, the task outputs three intermediate TSV files (features, barcodes, and matrix) used for downstream empty droplet detection and Loom matrix generation. 


#### 5. Calculate gene metrics

The [CalculateGeneMetrics](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) task uses [warp-tools](https://github.com/broadinstitute/warp-tools) to calculate summary metrics that help assess the quality of the data output each time this pipeline is run. These metrics are included in the output Loom matrix. A detailed list of these metrics is found in the [Optimus Count Matrix Overview](./Loom_schema.md).

#### 6. Calculate cell metrics

The [CalculateCellMetrics](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Metrics.wdl) task uses [warp-tools](https://github.com/broadinstitute/warp-tools) to calculate summary metrics that help assess the per-cell quality of the data output each time this pipeline is run. These metrics are included in the output Loom matrix. A detailed list of these metrics is found in the [Optimus Count Matrix Overview](./Loom_schema.md).


#### 7. Run emptyDrops

Empty droplets are lipid droplets that did not encapsulate a cell during 10x sequencing, but instead acquired cell-free RNA (secreted RNA or RNA released during cell lysis) from the solution in which the cells resided ([Lun, et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/?term=30902100). This ambient RNA can serve as a substrate for reverse transcription, leading to a small number of background reads. The Optimus pipeline calls the [RunEmptyDrops](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/RunEmptyDrops.wdl) task which uses the [DropletUtils](http://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package to flag CBs that represent empty droplets rather than cells. A cell will be flagged if it contains fewer than 100 molecules.
emptyDrops metrics for single-cell data (not single-nucleus; see note below) are stored in the columns of the output Loom matrix.  Read full details for all the metrics in the [Optimus Matrix Overview](./Loom_schema.md).

EmptyDrops relies on a visual knee point inflection (described in [Lun et al. (2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)) to differentiate ambient-like cells from empty droplets. **If the single-cell data (counting mode set to `sc_rna`) does not produce a knee point inflection when running emptyDrops, the Loom columns for emptyDrops data will contain "NA"s.**

:::warning RunEmptyDrops output not included for single-nucleus data
Often snRNAseq data does not produce a visual knee point inflection when running emptyDrops and the tool can not accurately distinguish ambient-like cells from empty droplets. For this reason, emptyDrops is not used if Optimus counting_mode is set to `sn_rna`, and the output Loom matrix will not contain emptyDrops metrics.
:::

#### 8.  Matrix construction

The [OptimusLoomGeneration](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/LoomUtils.wdl) task uses a custom python script to merge the converted STARsolo count matrix, the emptyDrops results, and the cell and gene metrics into a Loom-formatted cell-by-gene matrix. **These counts are raw and unfiltered.**

Read full details for all the metrics in the [Optimus Count Matrix Overview](./Loom_schema.md).

The type of gene counts in the Loom will vary depending on the Optimus workflow counting_mode. If running single-cell data (sc_rna mode), the counts will include only exonic gene counts. 

If running single-nucleus data (sn_rna mode), the counts in the main matrix will be whole transcript and, when `count_exons` is set to true, the counts in an additional layer will be exonic. Using the `count_exons` parameter will cause the Loom matrix to have additional columns (cell barcodes) due to the difference in STARsolo counting mode.

You can determine which type of counts are in the loom by looking at the global attribute `expression_data_type`.

For sn_rna mode, you can also access whole transcript and exonic counts using Loompy's `layers()` method. For example, `loompy.connect.layers[“”]` will return the whole transcript counts from the output Loom file. Similarly, `loompy.connect.layers[“exon_counts”]` will return the exonic counts from the output Loom. 


#### 9. Outputs

Output files of the pipeline include:

1. Cell x Gene unnormalized, but UMI-corrected, count matrices in Loom format.
2. Unfiltered, sorted BAM file with barcode and downstream analysis [tags](./Bam_tags.md).
3. Cell metadata, including cell metrics.
4. Gene metadata, including gene metrics.

The following table lists the output files produced from the pipeline. For samples that have sequenced over multiple lanes, the pipeline will output one merged version of each listed file.

| Output Variable Name | Filename, if applicable | Output Type |Output Format |
| ------ |------ | ------ | ------ |
| pipeline_version_out | N/A | Version of the processing pipeline run on this data. | String |
| bam | `<input_id>.bam` | Aligned BAM | BAM |
| matrix | `<input_id>_sparse_counts.npz` | Converted sparse matrix file from the MergeStarOutputs task. | NPZ |
| matrix_row_index | `<input_id>_sparse_counts_row_index.npy` | Index of cells in count matrix. | NPY |
| matrix_col_index | `<input_id>_sparse_counts_col_index.npy` | Index of genes in count matrix. | NPY |
| cell_metrics | `<input_id>.cell-metrics.csv.gz` | Cell metrics | Compressed CSV | Matrix of metrics by cells. |
| gene_metrics | `<input_id>.gene-metrics.csv.gz` | Gene metrics | Compressed CSV | Matrix of metrics by genes. |
| cell_calls | empty_drops_result.csv | emptyDrops results from the RunEmptyDrops task. | CSV |
| loom_output_file | `<input_id>.loom` | Loom | Loom | Loom file with count data (exonic or whole transcript depending on the counting_mode) and metadata. | N/A |

The Loom matrix is the default output. See the [create_loom_optimus.py](https://github.com/broadinstitute/warp/blob/master/dockers/skylab/pytools/tools/create_loom_optimus.py) for the detailed code. This matrix contains the unnormalized (unfiltered), UMI-corrected count matrices, as well as the gene and cell metrics detailed in the [Optimus Count Matrix Overview](./Loom_schema.md).

#### Try the Optimus matrix with community tools
The matrix is compatible with multiple downstream community analysis tools, including [Seurat](https://satijalab.org/seurat/index.html), [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Cumulus](https://cumulus.readthedocs.io/en/latest/index.html), and [Pegasus](https://pegasus.readthedocs.io/en/stable/#). To try a tutorial using the Optimus matrix with these tools, register for the open-source platform [Terra](https://app.terra.bio) and then navigate to the public [Intro-to-HCA-data-on-Terra workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/Intro-to-HCA-data-on-Terra). You can also view the accompanying [step-by-step guide](https://support.terra.bio/hc/en-us/articles/360060041772) without registration.


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


## Citing the Optimus pipeline
Please identify the pipeline in your methods section using the Optimus Pipeline's [SciCrunch resource identifier](https://scicrunch.org/scicrunch/Resources/record/nlx_144509-1/SCR_018908/resolver?q=SCR_018908&l=SCR_018908).
* Ex: *Optimus Pipeline (RRID:SCR_018908)*

## Consortia support
This pipeline is supported and used by the [Human Cell Atlas](https://www.humancellatlas.org/) (HCA) project and the [BRAIN Initiative Cell Census Network](https://biccn.org/) (BICCN). 

Each consortium may use slightly different reference files for data analysis or have different post-processing steps. Learn more by reading the [Consortia Processing](./consortia-processing.md) overview.

If your organization also uses this pipeline, we would like to list you! Please reach out to us by contacting [the WARP team](mailto:warp-pipelines-help@broadinstitute.org).

## Feedback

Please help us make our tools better by contacting [the WARP team](mailto:warp-pipelines-help@broadinstitute.org) for pipeline-related suggestions or questions.


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
For Read 1 sequences, the only minimum requirement is that reads are the combined lengths of the CB and UMIs (which will vary between 10x V1, V2, and V3 chemistry).

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


