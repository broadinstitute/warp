version 1.0

#WORKFLOW DEFINITION   
workflow VUMCBamToCram {
  input {
    File input_bam
    File ref_fasta
    String sample_name
  }

  call BamToCram {
    input:
      input_bam = input_bam,
      ref_fasta = ref_fasta,
      sample_name = sample_name,
  }

  output {
    File output_cram = BamToCram.output_cram
    File output_cram_index = BamToCram.output_cram_index
    File output_cram_md5 = BamToCram.output_cram_md5
  }
}

# TASK DEFINITIONS
task BamToCram {
  input {
    File input_bam
    File ref_fasta
    String sample_name
    Int machine_mem_gb = 4
    Float additional_disk_size = 20
    Float disk_multiplier = 2.5
  }
  Int disk_size = ceil(size(input_bam, "GB")) * disk_multiplier + additional_disk_size

  String output_name = "~{sample_name}.cram"

  command <<<
    samtools view -T ~{ref_fasta} -C -o ~{output_name} ~{input_bam}
    samtools index ~{output_name}
    md5sum ~{output_name} > ~{output_name}.md5
  >>>

  runtime {
    docker: "staphb/samtools:1.17"
    preemptible: 3
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_cram = "~{output_name}"
    File output_cram_index = "~{output_name}.crai"
    File output_cram_md5 = "~{output_name}.md5"
  }
}
