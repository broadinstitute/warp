version 1.0


import "../tasks/broad/DragmapAlignment.wdl" as DragmapAlignment
import "../structs/dna_seq/DNASeqStructs.wdl"
import "../verification/VerifyTasks.wdl"


workflow DragmapAlignCompare {

  input {
    Array[File] input_bams_1
    Array[File] input_bams_2
  }


  scatter (idx in range(length(input_bams_1))) {
    File truth_bam = input_bams_1[idx]
    File test_bam = input_bams_2[idx]

    call VerifyTasks.CompareBams {
      input:
        test_bam = test_bam,
        truth_bam =test_bam,
        lenient_header = false
    }
  }

  output {
  }
  meta {
    allowNestedInputs: true
  }
}
