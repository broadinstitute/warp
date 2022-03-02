version 1.0


import "../tasks/broad/DragmapAlignment.wdl" as DragmapAlignment
import "../structs/dna_seq/DNASeqStructs.wdl"
import "../verification/VerifyTasks.wdl"


workflow DragmapAlignCompare {

  input {
    Array[File] test_bams
    Array[File] truth_bams
  }


  scatter (idx in range(length(truth_bams))) {
    File truth_bam = truth_bams[idx]
    File test_bam = test_bams[idx]

    call VerifyTasks.CompareLargeBamFiles {
      input:
        test_bam = test_bam,
        truth_bam = truth_bam,
    }
  }

  output {
  }
  meta {
    allowNestedInputs: true
  }
}
