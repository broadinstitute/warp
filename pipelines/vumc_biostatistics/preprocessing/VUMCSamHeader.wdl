version 1.0

#WORKFLOW DEFINITION   
workflow VUMCSamHeader {
  input {
    File input_cram
    String sample_name
  }

  call SamHeader {
    input:
      input_cram = input_cram,
      sample_name = sample_name,
  }

  output {
    File sam_header_file = SamHeader.sam_header_file
  }
}

# TASK DEFINITIONS
# Validate a cram using Picard ValidateSamFile
task SamHeader {
  input {
    String input_cram
    String sample_name
  }

  String output_name = "~{sample_name}.header.txt"

  command <<<
    samtools view -H ~{input_cram} > ~{output_name}
  >>>

  runtime {
    docker: "staphb/samtools:1.17"
    preemptible: 3
    memory: "10 GB"
    disks: "local-disk 5 HDD"
  }

  output {
    File sam_header_file = "~{output_name}"
  }
}
