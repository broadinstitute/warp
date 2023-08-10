version 1.0

#WORKFLOW DEFINITION
workflow VUMCCramReads {
  input {
    File input_cram
    String sample_name
    String samtools_docker = "staphb/samtools:latest"
  }
  
  call CountCRAM {
    input:
      input_cram = input_cram,
      docker = samtools_docker,
  }

  output {
    Int unmapped_reads = CountCRAM.unmapped_reads
    Int mapped_reads = CountCRAM.mapped_reads
  }
}

# TASK DEFINITIONS
# Counting the number of mapped and unmapped reads in CRAM file using samtools
task CountCRAM{
  input{
    # Command parameters
    File input_cram

    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }

  Int disk_size = ceil(size(input_cram, "GB")) + addtional_disk_space_gb
  String NumUnmapped = "Unmapped.txt"
  String NumMapped = "Mapped.txt"

  command <<<
    # new samtools will output primary mapped line, we need to ignore it
    samtools flagstat ~{input_cram} | grep "mapped (" | grep -v "primary" | cut -f1 -d' ' > ~{NumMapped}

    # We don't need reference file for counting unmapped reads
    samtools view -c -f4 ~{input_cram} > ~{NumUnmapped}
  >>>

  runtime{
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output{
   Int unmapped_reads = read_int("~{NumUnmapped}")
   Int mapped_reads = read_int("~{NumMapped}")
  }
}


