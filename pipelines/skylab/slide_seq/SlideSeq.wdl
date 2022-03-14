version 1.0

import "../../../tasks/skylab/StarAlign.wdl" as StarAlign
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
        String read_structure
        String tar_star_reference
        String whitelist
        String output_bam_basename
    }

    parameter_meta {
        r1_fastq: "Array of Read 1 FASTQ files - forward read, contains cell barcodes and molecule barcodes"
        r2_fastq: "Array of Read 2 FASTQ files - reverse read, contains cDNA fragment generated from captured mRNA"
        i1_fastq: "(optional) Array of i1 FASTQ files - index read, for demultiplexing of multiple samples on one flow cell."
        sample_id: "Name of sample matching this file, inserted into read group header"
        tar_star_reference: ""
        whitelist: ""
        read_structure: ""
        output_bam_basename: ""
    }

    call FastqProcessing.FastqProcessingSlideSeq as SplitFastq {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            i1_fastq = i1_fastq,
            read_structure = read_structure,
            sample_id = sample_id
    }

    call StarAlign.STARsoloFastqSlideSeq as STARsoloFastqSlideSeq {
        input:
            r1_fastq = SplitFastq.fastq_R1_output,
            r2_fastq = SplitFastq.fastq_R2_output,
            white_list = whitelist,
            tar_star_reference = tar_star_reference,
            output_bam_basename = output_bam_basename
    }

    output {
        File bam_output = STARsoloFastqSlideSeq.bam_output
        File alignment_log = STARsoloFastqSlideSeq.alignment_log
        File general_log = STARsoloFastqSlideSeq.general_log
        File barcodes = STARsoloFastqSlideSeq.barcodes
        File features = STARsoloFastqSlideSeq.features
        File matric = STARsoloFastqSlideSeq.matrix

    }
}