version 1.0

import "../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineTasks.wdl" as Tasks
import "../../tasks/wdl/Alignment.wdl" as AlignmentTasks
import "../../structs/dna_seq/UltimaGenomicsWholeGenomeGermlineStructs.wdl" as Structs

workflow AlignmentAndMarkDuplicates {
  input {
    AlignmentReferences alignment_references
    References references

    Boolean is_cram
    Array[File] input_cram_bam
    Float rsq_threshold
    Int reads_per_split
    String base_file_name_sub
    Boolean save_bam_file
  }

  Int compression_level = 2

  if (is_cram) {
    scatter(input_cram in input_cram_bam) {
      call Tasks.SplitCram as SplitInputCram{
        input:
          input_cram_bam = input_cram,
          base_file_name = base_file_name_sub,
          reads_per_file = reads_per_split
      }
    }
  }
  if (!is_cram) {
    scatter(input_bam in input_cram_bam) {
      call AlignmentTasks.SamSplitter as SplitInputBam {
        input:
          input_bam         = input_bam,
          compression_level = compression_level,
          n_reads           = reads_per_split
      }
    }
  }

  Array[File] split_outputs = flatten(select_first([SplitInputCram.split_outputs, SplitInputBam.split_bams, []]))

  scatter(split_chunk in split_outputs) {
    Float split_chunk_size = size(split_chunk, "GB")

    call Tasks.ConvertCramOrBamToUBam as ConvertToUbam {
      input:
        input_file        = split_chunk,
        split_chunk_size  = split_chunk_size,
        base_file_name    = base_file_name_sub
    }
    
    call Tasks.SamToFastqAndBwaMemAndMba{
      input :
        input_bam            = ConvertToUbam.unmapped_bam,
        alignment_references = alignment_references,
        references           = references
    }
  }
  # Add one local ssd disk (375 GB) when spare disk is < 50 GB
  # the actual size that will be required to google will be a multiple of 375
  Int md_disk_size = 4 * ceil(size(SamToFastqAndBwaMemAndMba.output_bam, "GB"))
  Int rounded_disk_size = ceil(md_disk_size/375) * 375
  Int total_md_disk_size = if (ceil(md_disk_size/375) * 375 - 4 * ceil(size(SamToFastqAndBwaMemAndMba.output_bam, "GB"))) < 50 then 4 * ceil(size(SamToFastqAndBwaMemAndMba.output_bam, "GB")) + 375 else 4 * ceil(size(SamToFastqAndBwaMemAndMba.output_bam, "GB"))
  Int mapped_bam_size_local_ssd = if total_md_disk_size < 3000 then total_md_disk_size else 9000

  call Tasks.MarkDuplicatesSpark {
    input:
      input_bams          = SamToFastqAndBwaMemAndMba.output_bam,
      output_bam_basename = base_file_name_sub + ".aligned.sorted.duplicates_marked",
      save_bam_file       = save_bam_file,
      disk_size_gb        = mapped_bam_size_local_ssd
  }

  output {
    File output_bam = MarkDuplicatesSpark.output_bam
    File output_bam_index = MarkDuplicatesSpark.output_bam_index
    File? optional_output_bam = MarkDuplicatesSpark.optional_output_bam
    File? optional_output_bam_index = MarkDuplicatesSpark.optional_output_bam_index
  }
}