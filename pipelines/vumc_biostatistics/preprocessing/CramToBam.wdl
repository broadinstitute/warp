version 1.0

#WORKFLOW DEFINITION   
workflow VUMCBamToCram {
  input {
    File input_cram
    File? input_cram_index

    File? old_ref_fasta
    File? old_ref_fasta_index

    File ref_fasta
    File ref_fasta_index

    String sample_name
  }

  call BamToCram {
    input:
      input_cram = input_cram,
      input_cram_index = input_cram_index,
      old_ref_fasta = old_ref_fasta,
      old_ref_fasta_index = old_ref_fasta_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sample_name = sample_name,
  }

  output {
    File output_bam = CramToBam.output_bam
    File output_bam_index = CramToBam.output_bam_index
    File output_bam_md5 = CramToBam.output_bam_md5
  }
}

# TASK DEFINITIONS
task BamToCram {
  input {
    File input_cram
    File? input_cram_index

    File? old_ref_fasta
    File? old_ref_fasta_index

    File ref_fasta
    File ref_fasta_index

    String sample_name

    Int machine_mem_gb = 4
    Int additional_disk_size = 20
    Float disk_multiplier = 2.5
  }
  Int disk_size = ceil(size(input_cram, "GB") * disk_multiplier + additional_disk_size)

  String output_name = "~{sample_name}.bam"

  command <<<
    samtools view -b -T ~{ref_fasta} -o ~{output_name} ~{input_cram}
    samtools index ~{output_name}
    md5sum ~{output_name} > ~{output_name}.md5
  >>>