version 1.0

import "https://raw.githubusercontent.com/broadinstitute/warp/np_add_starsolo_slideseq/tasks/skylab/StarAlign.wdl" as StarAlign

workflow SlideSeq {

  input {

    # Sequencing data inputs
    File ubam

    # organism reference parameters
    File tar_star_reference
    File annotations_gtf

    # 10x parameters
    File whitelist
    String counting_mode = "sc_rna"
    String output_bam_basename

  }

call StarAlign.STARsoloFastq as STARsoloFastqSlideSeq {
    input:
      ubam = ubam,
      white_list = whitelist,
      tar_star_reference = tar_star_reference,
      counting_mode = counting_mode,
      output_bam_basename = output_bam_basename
      }

 output {

    File bam = STARsoloFastqSlideSeq.bam_output

  }
}