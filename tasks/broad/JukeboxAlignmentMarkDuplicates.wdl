version 1.0

import "JukeboxTasks.wdl" as Tasks
import "Alignment.wdl" as AlignmentTasks
import "../../structs/dna_seq/JukeboxStructs.wdl" as Structs

workflow AlignmentAndMarkDuplicates {
  input {
    SampleInputs sample_inputs
    String base_file_name_sub
    Int? reads_per_split
    String crammer_docker
    Int preemptibles
    Int compression_level
    Int additional_disk
    String gitc_docker
    String gitc_path
    String gatk_markduplicates_docker
    String jukebox_vc_docker

    Boolean no_address
    Boolean parallel_no_address
    String dummy_input_for_call_caching
    
    Float rsq_threshold

    String bwa_commandline
    AlignmentReferences alignment_references

    String? mark_duplicates_extra_args

    Float bwa_disk_multiplier
    Int local_ssd_size
    Float ref_size

    File monitoring_script
  }

  # Get the size of the standard reference files as well as the additional reference files needed for BWA
  Float bwa_ref_size = ref_size + size(alignment_references.ref_alt, "GB") + size(alignment_references.ref_amb, "GB") + size(alignment_references.ref_ann, "GB") + size(alignment_references.ref_bwt, "GB") + size(alignment_references.ref_pac, "GB") + size(alignment_references.ref_sa, "GB")

  String detect_input_ending_from_file = if sub(sample_inputs.input_cram_bam_list[0], ".*\\.cram$", "is_cram") == "is_cram" then "is_cram" else "is_bam"
  String detect_input_ending = select_first([sample_inputs.override_input_ending, detect_input_ending_from_file])

  if (detect_input_ending == "is_cram") {

    scatter(input_cram in sample_inputs.input_cram_bam_list) {

      Int reads_per_cram = select_first([reads_per_split, 20000000])

      call Tasks.SplitCram as SplitInputCram{
        input:
          monitoring_script = monitoring_script,
          input_cram_bam = input_cram,
          base_file_name = base_file_name_sub,
          reads_per_file = reads_per_cram,
          docker = crammer_docker,
          preemptible_tries = preemptibles,
          no_address = false
      }
    }
    Array[Array[File]]? split_cram = SplitInputCram.split_outputs
  }

  if ( detect_input_ending == "is_bam" ) {

    scatter(input_bam in sample_inputs.input_cram_bam_list){

      File unmapped_bam = input_bam
      Float unmapped_bam_size = size(unmapped_bam, "GB")

      # Split bam into multiple smaller bams,
      # map reads to reference and recombine into one bam
      Int reads_per_file = 20000000
      call AlignmentTasks.SamSplitter as SplitInputBam {
        input:
          input_bam = unmapped_bam,
          n_reads = reads_per_file,
          preemptible_tries = preemptibles,
          compression_level = compression_level
      }

    }
    Array[Array[File]]? split_bam = SplitInputBam.split_bams
  }

  Array[File] split_outputs = flatten(select_first([split_cram, split_bam, []]))

  call Tasks.GetBwaVersion {
    input:
      docker = gitc_docker,
      gitc_path = gitc_path,
      no_address = no_address,
      dummy_input_for_call_caching = dummy_input_for_call_caching,
  }

  scatter(split_chunk in split_outputs) {

    Float split_chunk_size = size(split_chunk, "GB")

    String split_chunk_name = if sub(split_chunk, ".*\\.cram$", "cram") == "cram" then basename(split_chunk, ".cram") else  basename(split_chunk, ".bam")

    call Tasks.ConvertCramOrBamToUBam as ConvertToUbam {
      input:
        monitoring_script = monitoring_script,
        input_file = split_chunk,
        split_chunk_size = split_chunk_size,
        additional_disk = additional_disk,
        base_file_name = base_file_name_sub,
        preemptible_tries = preemptibles,
        docker = gatk_markduplicates_docker,
        no_address = parallel_no_address
    }

    call Tasks.FilterByRsq {
      input:
        monitoring_script = monitoring_script, 
        docker = jukebox_vc_docker, 
        input_bam = ConvertToUbam.unmapped_bam, 
        rsq_threshold = rsq_threshold,
        no_address = no_address, 
        preemptible_tries = preemptibles,
        additional_disk = additional_disk
    }

    File unmapped_bam = FilterByRsq.output_bam
    Float unmapped_bam_size = size(unmapped_bam, "GB")
    String current_name = basename(unmapped_bam, ".bam")

    call Tasks.SamToFastqAndBwaMemAndMba{
      input :
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = current_name,
        alignment_references = alignment_references,
        bwa_version = GetBwaVersion.version,
        compression_level = compression_level,
        preemptible_tries = preemptibles,
        docker = gitc_docker,
        monitoring_script=monitoring_script,
      # The merged bam can be bigger than only the aligned bam,
      # so account for the output size by multiplying the input size by 2.75.
        disk_size = ceil(unmapped_bam_size + bwa_ref_size + (bwa_disk_multiplier * unmapped_bam_size) + additional_disk),
        no_address = parallel_no_address
    }
  }

  Int bams_sum = ceil(size(SamToFastqAndBwaMemAndMba.output_bam, "GB"))
  Int md_disk_size = 4 * ceil(bams_sum)
  # add one local ssd disk (375 GB) when spare disk is < 50 GB
  # the actual size that will be required to google will be a multiple of 375
  Int rounded_disk_size = ceil(md_disk_size/local_ssd_size) * local_ssd_size
  Int total_md_disk_size = if (rounded_disk_size - md_disk_size) < 50 then md_disk_size + local_ssd_size else md_disk_size

  Int local_ssd_threshold = 3000
  Int mapped_bam_size_local_ssd = if total_md_disk_size < local_ssd_threshold then total_md_disk_size else 9000

  # increase memory size for large inputs
  Int increased_ram_size = 300
  Int mark_duplicates_memory_threshold = 600
  Int mark_duplicates_ram = if md_disk_size / 4 > mark_duplicates_memory_threshold then increased_ram_size else 208
  Int mark_duplicates_cpus = 32

  call Tasks.MarkDuplicatesSpark {
    input:
      input_bams = SamToFastqAndBwaMemAndMba.output_bam,
      output_bam_basename = base_file_name_sub + ".aligned.sorted.duplicates_marked",
      monitoring_script = monitoring_script,
      memory_gb = mark_duplicates_ram,
      disk_size_gb = mapped_bam_size_local_ssd,
      docker = gatk_markduplicates_docker,
      gitc_path = gitc_path,
      no_address = false,
      cpu = mark_duplicates_cpus,
      args = mark_duplicates_extra_args
  }

  output {
    File output_bam = MarkDuplicatesSpark.output_bam
    File output_bam_index = MarkDuplicatesSpark.output_bam_index
  }
}