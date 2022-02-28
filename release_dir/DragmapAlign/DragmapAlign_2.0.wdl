version 1.0


import "DragmapAlignment.wdl" as DragmapAlignment
import "DNASeqStructs.wdl"
import "VerifyTasks.wdl"


workflow DragmapAlign {

  input {
    Array[File] input_bams

    ReferenceFasta reference_fasta
    DragmapReference dragmap_reference

    Int compression_level
    Int preemptible_tries
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true

    String docker = "us.gcr.io/broad-gotc-prod/dragmap:1.1.0-1.2.1-2.26.4-1.11-1638901436"
    Int cpu = 16
    Float disk_multiplier = 8
    Int memory_mb = 40960
  }
  scatter (x in range(20)) {
    scatter (input_bam in input_bams) {
      String unmapped_bam_basename = basename(input_bam, ".unmapped.bam")
      call DragmapAlignment.SamToFastqAndDragmapAndMba as DragmapAlignFirstRun {
        input:
          input_bam = input_bam,
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
          reference_fasta = reference_fasta,
          dragmap_reference = dragmap_reference,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries,
          hard_clip_reads = hard_clip_reads,
          unmap_contaminant_reads = unmap_contaminant_reads,
          cpu_platform = "Intel Haswell"
      }

      call DragmapAlignment.SamToFastqAndDragmapAndMba as DragmapAlignSecondRun {
        input:
          input_bam = input_bam,
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
          reference_fasta = reference_fasta,
          dragmap_reference = dragmap_reference,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries,
          hard_clip_reads = hard_clip_reads,
          unmap_contaminant_reads = unmap_contaminant_reads,
          cpu_platform = "Intel Skylake"
      }

      call VerifyTasks.CompareLargeBamFiles {
        input:
          test_bam = DragmapAlignSecondRun.output_bam,
          truth_bam = DragmapAlignFirstRun.output_bam,
      }
    }
  }
  output {
    Array[Array[File]] first_bam = DragmapAlignFirstRun.output_bam
    Array[Array[File]] first_stderr_log = DragmapAlignFirstRun.dragmap_stderr_log
    Array[Array[File]] second_bam = DragmapAlignSecondRun.output_bam
    Array[Array[File]] second_stderr_log = DragmapAlignSecondRun.dragmap_stderr_log
  }
  meta {
    allowNestedInputs: true
  }
}
