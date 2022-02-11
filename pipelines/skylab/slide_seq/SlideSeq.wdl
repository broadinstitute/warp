version 1.0

import "https://raw.githubusercontent.com/broadinstitute/warp/np_add_starsolo_slideseq/tasks/skylab/StarAlign.wdl" as StarAlign

workflow SlideSeq {

  input {

    # Sequencing data inputs
    Array[File] r1_fastq
    Array[File] r2_fastq

    # organism reference parameters
    File tar_star_reference
    File annotations_gtf

    # 10x parameters
    File whitelist
    String output_bam_basename

    Int umi_length
    Int cell_barcode_length

  }

call StarAlign.STARsoloFastqSlideSeq as STARsoloFastqSlideSeq {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      white_list = whitelist,
      tar_star_reference = tar_star_reference,
      output_bam_basename = output_bam_basename
      }

 output {
    File bam = STARsoloFastqSlideSeq.bam_output
  }
}