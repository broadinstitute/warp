---
sidebar_position: 2
---

# Optimus v7.8.1 Methods

Below we provide an example methods section for a publication, separated into single-cell or single-nucleus use cases. For the complete pipeline documentation, see the [Optimus Overview](./README.md).

# Methods

## Single-cell (sc_rna mode)
Data preprocessing and count matrix construction were performed using the Optimus v7.8.1 pipeline (RRID:SCR_018908). Briefly, FASTQ files were partitioned by barcodes using fastqprocess. The files were then trimmed, aligned, UMI-corrected against the 10x Genomics barcodes whitelist, and converted to a raw count matrix using STARsolo v2.7.11a. CB correction was performed using the `--soloCBmatchWLtype 1MM_multi` parameter.

Reads were trimmed using the solo parameter `--clipAdapterType CellRanger4` and `--outFilterScoreMin 30` which matches read trimming performed by CellRanger4. Reads were then aligned to GENCODE mouse (M32) or human (V43) references in stranded mode. Genes were annotated using the STARsolo "Gene" COUNTING_MODE and UMIs were corrected with the `--soloUMIdedup 1MM_CR` parameter, which uses Cell Ranger's correction method. The resulting BAM was then used for cell and gene metric correction using the warp-tools TagSort tool. The STAR TSV outputs for gene counts, features, and barcodes were converted to numpy arrays for downstream empty droplet detection using DropletUtils v1.2.1 emptyDrops with the parameters `--fdr-cutoff 0.01 --emptydrops-niters 10000 --min-molecules 100 --emptydrops-lower 100.`


All cell and gene metrics (alignment, mitochondrial, and other QC metrics), count matrices, and emptyDrops results were aggregated into a final h5ad-formatted cell-by-gene matrix. The final outputs included the unfiltered h5ad and unfiltered (but tagged) BAM file.

An example of the pipeline and outputs is available on the [Terra HCA Optimus Pipeline Featured Workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline), and an additional pipeline overview is available in [WARP documentation.](https://broadinstitute.github.io/warp/docs/Pipelines/Optimus_Pipeline/README) Examples of genomic references, whitelists, and other inputs are available in the WARP repository [(see the example inputs).](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/optimus/example_inputs).

## Single-nucleus (sn_rna mode)

Data preprocessing and count matrix construction were performed using the Optimus v7.8.1 pipeline (RRID:SCR_018908). Briefly, FASTQ files were partitioned by barcodes using fastqprocess. The files were then trimmed, aligned, UMI-corrected against the 10x Genomics barcodes whitelist, and converted to a raw count matrix using STARsolo v2.7.11a. CB correction was performed using the `--soloCBmatchWLtype 1MM_multi` parameter.

Reads were trimmed using the solo parameter `--clipAdapterType CellRanger4` and `--outFilterScoreMin 30` which matches read trimming performed by CellRanger4. Reads were then aligned to GENCODE mouse (M32) or human (V43) references in stranded mode. Genes were annotated using the STAR "GeneFull_Ex50pAS" COUNTING_MODE and UMIs were corrected with the `--soloUMIdedup 1MM_CR`, which uses a Cell Ranger's correction method. The resulting BAM was then used for cell and gene metric correction using the warp-tools TagSort tool. The STAR TSV outputs for gene counts, features, and barcodes were converted to numpy arrays for downstream h5ad conversion. All cell and gene metrics (alignment, mitochondrial, and other QC metrics) and count matrices were aggregated into a final h5ad-formatted cell-by-gene matrix. The final outputs included the unfiltered h5ad and unfiltered (but tagged) BAM file.

An example of the pipeline and outputs is available on the [Terra HCA Optimus Pipeline Featured Workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline), and an additional pipeline overview is available in [WARP documentation](https://broadinstitute.github.io/warp/docs/Pipelines/Optimus_Pipeline/README). Examples of genomic references, whitelists, and other inputs are available in the WARP repository (see the [example inputs.](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/optimus/example_inputs)).


