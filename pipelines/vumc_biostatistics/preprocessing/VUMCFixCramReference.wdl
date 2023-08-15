version 1.0

#WORKFLOW DEFINITION   
workflow VUMCFixCramReference {
  input {
    File input_bam
    File? input_bam_index

    File? old_ref_fasta
    File? old_ref_fasta_index

    File ref_fasta
    File ref_fasta_index

    String sample_name
  }

  call FixReference {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      old_ref_fasta = old_ref_fasta,
      old_ref_fasta_index = old_ref_fasta_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sample_name = sample_name,
  }

  output {
    File output_cram = FixReference.output_cram
    File output_cram_index = FixReference.output_cram_index
    File output_cram_md5 = FixReference.output_cram_md5
  }
}

# TASK DEFINITIONS
task FixReference {
  input {
    File input_bam
    File? input_bam_index

    File? old_ref_fasta
    File? old_ref_fasta_index

    File ref_fasta
    File ref_fasta_index

    String sample_name

    Int machine_mem_gb = 4
    Int additional_disk_size = 20
    Float disk_multiplier = 2.5
  }
  Int disk_size = ceil(size(input_bam, "GB") * disk_multiplier + additional_disk_size)

  String output_name = "~{sample_name}.cram"

  command <<<

samtools view -h -T ~{ref_fasta} -C --no-PG --output-fmt-option embed_ref=2 ~{input_bam} | \
samtools view -h -T ~{ref_fasta} --no-PG - | grep -v '@PG' | \
samtools view -T ~{ref_fasta} -C --no-PG -o ~{output_name} -

samtools index ~{output_name}

md5sum ~{output_name} > ~{output_name}.md5.txt

  >>>

  runtime {
    docker: "evolbioinfo/samtools:v1.18"
    preemptible: 3
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_cram = "~{output_name}"
    File output_cram_index = "~{output_name}.crai"
    File output_cram_md5 = "~{output_name}.md5.txt"
  }
}
