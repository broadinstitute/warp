---
sidebar_position: 2
---

# scATAC v1.2.0 Methods

Below we provide a sample methods section for a publication. For the complete pipeline documentation, see the [scATAC Overview](./README.md).

## Methods

Preprocessing and cell by bin matrix construction were performed using the scATAC v1.2.0 Pipeline (RRID:SCR_018919).

Prior to processing, paired-end FASTQ files were modified with a custom python script (available at https://github.com/r3fang/SnapTools/blob/master/snaptools/dex_fastq.py) so that readnames were appended with cell barcodes.

The appended reads were then aligned to the Hg38 genomic reference using BWA v0.7.17. The resulting aligned BAM was converted into fragments and filtered using the [SnapTools v1.4.7](https://github.com/r3fang/SnapTools) SnapPre function with default parameters.

The snap-add-bmat function was then used to add cell-by-bin matrices to the resulting Snap file. 10 kb was selected as the default value for the bin size based on the SnapTools recommendation for mid-size datasets. It can be changed by specifying the desired size as an input to this workflow.

Custom python scripts were then used to make a [GA4GH-compliant BAM](https://github.com/broadinstitute/warp/blob/master/dockers/skylab/snaptools/makeCompliantBAM.py) and to export select Snap file metrics to [individual text files](https://github.com/broadinstitute/warp/blob/master/dockers/skylab/snap-breakout/breakoutSnap.py).

An example of the pipeline and its outputs is available on [Terra](https://app.terra.bio/#workspaces/brain-initiative-bcdc/SnapATAC_Pipeline) and more documentation can be found at [here](./README.md).
Examples of genomic reference files and other inputs can be found in the pipeline’s [example JSON](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/scATAC/example_inputs/human_example.json).
