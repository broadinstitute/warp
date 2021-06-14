---
sidebar_position: 2
---

# Optimus v4.2.4 Methods

Below we provide a sample methods sections for a publication, separated into single-cell or single-nucleus use cases. For the complete pipeline documentation, see the [Optimus Overview](./README.md).

# Methods

## Single-cell (sc_rna mode)
Data preprocessing and count matrix construction were performed using the Optimus v4.2.4 Pipeline (RRID:SCR_018908). Briefly, FASTQ files were converted to unaligned BAM (uBAM) using Picard v2.10.10 and reads were appended with raw UMI and corrected cell barcode sequences using Single Cell Tools (sctools) v0.3.12 and the 10x Genomics barcodes whitelist, allowing for up to one edit distance (Levenshtein distance).

uBAMs were then aligned to GENCODE mouse (M21) or human (V27) references using STAR v2.5.3a with default parameters in addition to

```
--BAM unsorted --outSAMattributes all --outSAMunmapped --readFilestype SAM SE
```

Genes were annotated and reads were tagged with Drop-seq Tools v1.12 using the TagReadwithGeneExon function.

UMIs were then corrected and duplicate reads marked using UMI-tools v0.0.1 with default parameters in addition to

```
--extract-umi-method=tag --umi-tag UR --cell-tag CB --gene-tag GE --umi-group-tag UB --per-gene --per-cell --no-sort-output
```

All reads (UMI-corrected, duplicate, and untagged) were merged into a single BAM file and tagged. Gene and cell-specific metrics were calculated using the sctools v0.3.11 functions `CalculateGeneMetrics` and `CalculateCellMetrics`.

Empty droplets were identified, but not removed to enable downstream filtering, using the DropletUtils v.1.2.1 with

```
--fdr-cutoff 0.01 --emptydrops-niters 10000 --min-molecules 100 --emptydrops-lower 100
```

UMI-aware count matrices for exon-only alignments were produced using the sctools v0.3.11.
All cell and gene metrics (alignment, mitochondrial, and other QC metrics), count matrices and DropletUtils results were then aggregated into a final Loom file for downstream processing. The final outputs included the unfiltered Loom and unfiltered (but tagged) BAM files.

An example of the pipeline and outputs is available on the Terra HCA Optimus Pipeline Featured Workspace (https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline), and additional documentation is available on GitHub (https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/README.md). Examples of genomic references, whitelists, and other inputs are available in the Skylab repository (see example JSONs).

## Single-nucleus (sn_rna mode)

Data preprocessing and count matrix construction were performed using the Optimus v4.1.7 Pipeline (RRID:SCR_018908). Briefly, FASTQ files were converted to unaligned BAM (uBAM) using Picard v2.10.10 and reads were appended with raw UMI and corrected cell barcode sequences using Single Cell Tools (sctools) v0.3.11 and the 10x Genomics barcodes whitelist, allowing for up to one edit distance (Levenshtein distance).

uBAMs were aligned to GENCODE mouse (M21) or human (V27) references using STAR v2.5.3a with default parameters in addition to

```
--BAM unsorted --outSAMattributes all --outSAMunmapped --readFilestype SAM SE
```

Genes were annotated and reads were tagged with Drop-seq Tools v2.3.0 using TagReadWithGeneFunction.

UMIs were then corrected using UMI-tools v0.0.1 with default parameters in addition to

```
--extract-umi-method=tag --umi-tag UR --cell-tag CB --gene-tag GE --umi-group-tag UB --per-gene --per-cell --no-sort-output
```

All reads (UMI-corrected, duplicate, and untagged) were merged into a single BAM file. Gene and cell-specific metrics were calculated using the sctools v0.3.11 functions `CalculateGeneMetrics` and `CalculateCellMetrics`. UMI-aware count matrices for all alignments (introns, exons, UTRs) were produced using the sctools v0.3.11.

All cell and gene metrics (alignment, mitochondrial, and other QC metrics), annotations, and count matrices were aggregated into a final Loom file for downstream processing. The final outputs included the unfiltered Loom and unfiltered (but tagged) BAM files.

An example of the pipeline and outputs is available on the [Terra HCA Optimus Pipeline Featured Workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline), and additional documentation is available on [GitHub](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/README.md). Examples of genomic references, whitelists, and other inputs are available in the Skylab repository (see the *_example.json files at [here](https://github.com/HumanCellAtlas/skylab/tree/master/pipelines/optimus).
