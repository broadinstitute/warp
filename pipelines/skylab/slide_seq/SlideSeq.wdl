version 1.0

import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing

## Copyright Broad Institute, 2021
##
## This WDL pipeline implements data processing for RNA with UMIs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.


workflow SlideSeq{

    String pipeline_version = "0.0.1"

    input {
        Array[File] r1_fastq
        Array[File] r2_fastq
        String sample_id

        Int? cell_barcode_length
        Int? umi_length
    }

    parameter_meta {
        r1_fastq: "Array of Read 1 FASTQ files"
        r2_fastq: "Array of Read 2 FASTQ files"
        sample_id: "Id for the sample being processed"
        cell_barcode_length: "Number of cell barcode base pairs in the Read 1 FASTQ"
        umi_length: "Number of UMI base pairs in the Read 1 FASTQ"
    }

}