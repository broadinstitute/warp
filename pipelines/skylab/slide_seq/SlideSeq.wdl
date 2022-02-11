version 1.0

import "../../../tasks/skylab/StarAlign.wdl" as StarAlign

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

        #TODO eventually get rid of these inputs and just use read_structure. But to do this, we need to fix fastqprocess first
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