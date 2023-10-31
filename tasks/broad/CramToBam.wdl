version 1.0

workflow CramToBam {

  input {
    File    input_cram
    File    ref_fasta
    File    ref_fasta_index
    String  output_basename
  }

  # Convert the final merged recalibrated BAM file to CRAM format
  call ConvertToBam {
    input:
      input_cram      = input_cram,
      ref_fasta       = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = output_basename
  }

  output {
     File output_bam = ConvertToBam.output_bam
     File output_bam_index = ConvertToBam.output_bam_index
  }

  meta {
    allowNestedInputs: true
  }
  
}

# Convert CRAM file to BAM format
task ConvertToBam {
  input {
    File input_cram
    File ref_fasta
    File ref_fasta_index
    String output_basename
  }

  command <<<
    set -e
    set -o pipefail

    samtools view -b -o ~{output_basename}.bam -T ~{ref_fasta} ~{input_cram}

    samtools index ~{output_basename}.bam
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    preemptible: 3
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk 200 HDD"
  }
  output {
    File output_bam = "~{output_basename}.bam"
    File output_bam_index = "~{output_basename}.bam.bai"
  }
}