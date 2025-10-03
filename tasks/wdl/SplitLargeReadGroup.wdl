version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements a split of large readgroups for human whole-genome and exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../tasks/wdl/Alignment.wdl" as Alignment
import "../../tasks/wdl/DragmapAlignment.wdl" as DragmapAlignment
import "../../tasks/wdl/BamProcessing.wdl" as Processing
import "../../tasks/wdl/Utilities.wdl" as Utils
import "../../structs/dna_seq/DNASeqStructs.wdl" as Structs

workflow SplitLargeReadGroup {

  input {
    File input_bam

    String bwa_commandline
    String output_bam_basename

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta
    DragmapReference? dragmap_reference

    Int compression_level
    Int preemptible_tries
    Int reads_per_file = 48000000
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
    Boolean use_bwa_mem = true
    Boolean allow_empty_ref_alt = false
  }

  call Alignment.SamSplitter as SamSplitter {
    input :
      input_bam = input_bam,
      n_reads = reads_per_file,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    Float current_unmapped_bam_size = size(unmapped_bam, "GiB")
    String current_name = basename(unmapped_bam, ".bam")

    if (use_bwa_mem) {
      call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          output_bam_basename = current_name,
          reference_fasta = reference_fasta,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries,
          hard_clip_reads = hard_clip_reads,
          unmap_contaminant_reads = unmap_contaminant_reads,
          allow_empty_ref_alt = allow_empty_ref_alt
      }
    }
    if (!use_bwa_mem) {
      call DragmapAlignment.SamToFastqAndDragmapAndMba as SamToFastqAndDragmapAndMba {
        input:
          input_bam = unmapped_bam,
          output_bam_basename = current_name,
          reference_fasta = reference_fasta,
          dragmap_reference = select_first([dragmap_reference]),
          compression_level = compression_level,
          preemptible_tries = preemptible_tries,
          hard_clip_reads = hard_clip_reads,
          unmap_contaminant_reads = unmap_contaminant_reads
      }
    }

    File output_bam = select_first([SamToFastqAndBwaMemAndMba.output_bam, SamToFastqAndDragmapAndMba.output_bam])
  }

  call Processing.GatherUnsortedBamFiles as GatherMonolithicBamFile {
    input:
      input_bams = output_bam,
      total_input_size = size(output_bam, "GiB"),
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }
  output {
    File aligned_bam = GatherMonolithicBamFile.output_bam
  }
  meta {
    allowNestedInputs: true
  }
}
