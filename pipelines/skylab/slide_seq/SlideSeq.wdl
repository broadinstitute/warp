version 1.0

import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing

## Copyright Broad Institute, 2022
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
        Array[File]? i1_fastq
        String sample_id

        Int cell_barcode_length
        Int umi_length
    }

    parameter_meta {
        r1_fastq: "Array of Read 1 FASTQ files - forward read, contains cell barcodes and molecule barcodes"
        r2_fastq: "Array of Read 2 FASTQ files - reverse read, contains cDNA fragment generated from captured mRNA"
        i1_fastq: "(optional) Array of i1 FASTQ files - index read, for demultiplexing of multiple samples on one flow cell."
        sample_id: "Name of sample matching this file, inserted into read group header"
        cell_barcode_length: "Number of cell barcode base pairs in the Read 1 FASTQ"
        umi_length: "Number of UMI base pairs in the Read 1 FASTQ"
    }

    call FastqProcessing.FastqProcessingSlidSeq as SplitFastq {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            i1_fastq = i1_fastq,
            umi_length = umi_length,
            cell_barcode_length = cell_barcode_length,
            sample_id = sample_id
    }

    output {

    }

}