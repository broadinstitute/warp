# ATAC v2.3.1 Methods

# Methods

Data preprocessing and analysis for 10x chromatin accessibility was performed using the ATAC workflow v2.3.1 (RRID:SCR_025042). Briefly, FASTQ files were processed with a custom tool fastqprocess which corrects cell barcodes against a reference whitelist and splits reads by barcode to enable processing parallelization. Adaptor sequences were then removed from reads using Cutadapt v4.4. Reads were then aligned to the reference genome using BWA-MEM2 v2.2.1 with default parameters, which outputs corrected barcodes to a BAM in the CB:Z tag. The resulting BAM was then processed with SnapATAC2 v2.7.0 to produce a fragment file, index, and h5ad containing fragments as well as per-barcode quality metrics.

An overview of the pipeline is available in [WARP Documentation](https://broadinstitute.github.io/warp/docs/Pipelines/ATAC/README) and examples of genomic references, whitelists, and other inputs are available in the [WARP repository](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/multiome/test_inputs).