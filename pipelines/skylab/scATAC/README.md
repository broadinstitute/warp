
## Announcing a new site for WARP documentation!

scATAC documentation has moved! Read more about the [scATAC Pipeline](https://broadinstitute.github.io/warp/docs/Pipelines/Single_Cell_ATAC_Seq_Pipeline/README) on the new [WARP documentation site](https://broadinstitute.github.io/warp/)!

### scATAC summary


The scATAC Pipeline was developed by the Broad DSP Pipelines team to process single cell/nucleus ATAC-seq datasets. The pipeline is based on the [SnapATAC pipeline](https://github.com/r3fang/SnapATAC) described by [Fang et al. (2019)](https://www.biorxiv.org/content/10.1101/615179v2.full). Overall, the pipeline uses the python module [SnapTools](https://github.com/r3fang/SnapTools) to align and process paired reads in the form of FASTQ files. It produces an hdf5-structured Snap file that includes a cell-by-bin count matrix. In addition to the Snap file, the final outputs include a GA4GH compliant aligned BAM and QC metrics.
