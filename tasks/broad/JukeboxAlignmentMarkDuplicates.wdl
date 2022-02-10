version 1.0

import "JukeboxTasks.wdl" as Tasks
import "Alignment.wdl" as AlignmentTasks
import "../../structs/dna_seq/JukeboxStructs.wdl" as Structs

workflow AlignmentAndMarkDuplicates {
  input {
    SampleInputs sample_inputs
    AlignmentReferences alignment_references

    Float ref_size
    Float rsq_threshold
    Int additional_disk
    String base_file_name_sub
    String dummy_input_for_call_caching

    String? mark_duplicates_extra_args
    Int? reads_per_split
  }

  Int compression_level = 2
  Float bwa_disk_multiplier = 2.5
  Float bwa_ref_size = ref_size
                      + size(alignment_references.ref_alt, "GB") 
                      + size(alignment_references.ref_amb, "GB") 
                      + size(alignment_references.ref_ann, "GB") 
                      + size(alignment_references.ref_bwt, "GB") 
                      + size(alignment_references.ref_pac, "GB") 
                      + size(alignment_references.ref_sa, "GB")
  String detect_input_ending_from_file = if sub(sample_inputs.input_cram_bam_list[0], ".*\\.cram$", "is_cram") == "is_cram" then "is_cram" else "is_bam"
  String detect_input_ending = select_first([sample_inputs.override_input_ending, detect_input_ending_from_file])

  # Values to dynamically determine disk & mem for MarkDuplicatesSpark
  Int local_ssd_size = 375
  Int local_ssd_threshold = 3000
  Int increased_ram_size = 300
  Int mark_duplicates_memory_threshold = 600
  Int mark_duplicates_cpus = 32

  if (detect_input_ending == "is_cram") {
    scatter(input_cram in sample_inputs.input_cram_bam_list) {
      call Tasks.SplitCram as SplitInputCram{
        input:
          input_cram_bam = input_cram,
          base_file_name = base_file_name_sub,
          reads_per_file = select_first([reads_per_split, 20000000])
      }
    }
  }

  if ( detect_input_ending == "is_bam" ) {
    scatter(input_bam in sample_inputs.input_cram_bam_list) {
      call AlignmentTasks.SamSplitter as SplitInputBam {
        input:
          input_bam         = input_bam,
          compression_level = compression_level,
          n_reads           = select_first([reads_per_split, 20000000])
      }
    }
  }

  Array[File] split_outputs = flatten(select_first([SplitInputCram.split_outputs, SplitInputBam.split_bams, []]))

  call Tasks.GetBwaVersion {
    input:
      dummy_input_for_call_caching = dummy_input_for_call_caching,
  }

  scatter(split_chunk in split_outputs) {
    Float split_chunk_size = size(split_chunk, "GB")
    
    call Tasks.ConvertCramOrBamToUBam as ConvertToUbam {
      input:
        input_file        = split_chunk,
        split_chunk_size  = split_chunk_size,
        additional_disk   = additional_disk,
        base_file_name    = base_file_name_sub,
    }

    call Tasks.FilterByRsq {
      input:
        input_bam       = ConvertToUbam.unmapped_bam, 
        rsq_threshold   = rsq_threshold,
        additional_disk = additional_disk
    }

    File unmapped_bam       = FilterByRsq.output_bam
    Float unmapped_bam_size = size(unmapped_bam, "GB")
    String current_name     = basename(unmapped_bam, ".bam")

    call Tasks.SamToFastqAndBwaMemAndMba{
      input :
        input_bam            = unmapped_bam,
        output_bam_basename  = current_name,
        alignment_references = alignment_references,
        bwa_version          = GetBwaVersion.version,
        disk_size            = ceil(unmapped_bam_size + bwa_ref_size + (bwa_disk_multiplier * unmapped_bam_size) + additional_disk),
    }
  }
  # Add one local ssd disk (375 GB) when spare disk is < 50 GB
  # the actual size that will be required to google will be a multiple of 375
  Int md_disk_size = 4 * ceil(size(SamToFastqAndBwaMemAndMba.output_bam, "GB"))
  Int rounded_disk_size = ceil(md_disk_size/local_ssd_size) * local_ssd_size
  Int total_md_disk_size = if (rounded_disk_size - md_disk_size) < 50 then md_disk_size + local_ssd_size else md_disk_size
  Int mapped_bam_size_local_ssd = if total_md_disk_size < local_ssd_threshold then total_md_disk_size else 9000
  Int mark_duplicates_ram = if md_disk_size / 4 > mark_duplicates_memory_threshold then increased_ram_size else 208


  call Tasks.MarkDuplicatesSpark {
    input:
      input_bams          = SamToFastqAndBwaMemAndMba.output_bam,
      output_bam_basename = base_file_name_sub + ".aligned.sorted.duplicates_marked",
      memory_gb           = mark_duplicates_ram,
      disk_size_gb        = mapped_bam_size_local_ssd,
      cpu                 = mark_duplicates_cpus,
      args                = mark_duplicates_extra_args
  }

  output {
    File output_bam = MarkDuplicatesSpark.output_bam
    File output_bam_index = MarkDuplicatesSpark.output_bam_index
  }
}