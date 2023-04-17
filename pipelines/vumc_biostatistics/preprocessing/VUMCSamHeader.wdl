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
    File input_cram
    String sample_name
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 10
  }
  Int disk_size = ceil(size(input_cram, "GB")) + addtional_disk_space_gb

  String output_name = "~{sample_name}.header.txt"

  command <<<
    samtools view -H ~{input_cram} > ~{output_name}
  >>>

  runtime {
    docker: "staphb/samtools:1.17"
    preemptible: 3
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File sam_header_file = "~{output_name}"
  }
}
