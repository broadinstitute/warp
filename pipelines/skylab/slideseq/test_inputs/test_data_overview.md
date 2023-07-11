# SlideSeq Test Data

The following lists the data used for workflow validation, including workflow engineering (plumbing) tests and scientific validation tests.

All unmodified FASTQ files are stored in a public Google bucket (gs://nemo-public/biccn-unbundled/grant/rf1_macosko/macosko/spatial_transcriptome/cellgroup/Slide-seq/mouse/raw/) and can be identified using their sample ID (ex. Puck_210817_11).

## Plumbing

* Puck_210817_11.mm10

The plumbing input FASTQ listed above was downsampled using a [custom python script](Plumbing/slideseq_downsample.py). First, the XY coordinates of a region of interest were selected by examining the total UMI counts per bead in a PDF image of the sample puck. A bead position TSV containing the XY coordinates for the puck's barcodes was used as input to the custom script and the TSV was cropped 2500 to 4000 in the x-dimension and 1000 to 4200 in the y-dimension. Next, 10 random cells were sampled for each 30 by 30 region in the cropped TSV. The barcodes for the 10 sampled cells were then used  as input to a second [custom script](https://github.com/broadinstitute/warp-tools/blob/master/fastqpreprocessing/src/samplefastq.cpp) from the [warp-tools GitHub repository](https://github.com/broadinstitute/warp-tools/tree/master). The output of this script was the new, downsampled FASTQ files.

## Scientific

* Puck_200501_23.mm10
* Puck_210727_12.mm10
* Puck_210727_17.mm10
* Puck_210817_07.mm10
* Puck_210817_11.mm10
* Puck_210817_12.mm10
* Puck_210817_13.mm10
* Puck_211013_03.mm10

The scientific inputs listed above were not downsampled.