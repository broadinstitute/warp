version 1.0

#WORKFLOW DEFINITION
workflow VUMCCramQC2 {
  input {
    Array[File] input_crams
    String sample_name
    String samtools_docker = "staphb/samtools:latest"
  }

  
 scatter (input_cram in input_crams){
  call CountCRAM {
    input:
      input_cram = input_cram,
      docker = samtools_docker,
  }
 }

 call SumUp {
  input: 
   sample_name = sample_name,
   mapped_files = CountCRAM.mapped_file,
   unmapped_files = CountCRAM.unmapped_file,
 }

output {
    Int unmapped_reads = SumUp.NumberUnmappedReads
    Int mapped_reads = SumUp.NumberMappedReads
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
   File unmapped_file = "~{NumUnmapped}"
   File mapped_file = "~{NumMapped}"
  }
  
}

task SumUp{
  input{
    # Command parameters
     String sample_name
     Array[File] unmapped_files
     Array[File] mapped_files
  }

  String FinalNumUnmapped = "${sample_name}_final_Unmapped.txt"
  String FinalNumMapped = "${sample_name}_final_Mapped.txt"

  command <<<
    # for array of file, we need to use sep to join them
    cat ~{sep=" " mapped_files} | awk '{ sum += $1 } END { print sum }' > ~{FinalNumMapped}
    cat ~{sep=" " unmapped_files} | awk '{ sum += $1 } END { print sum }' > ~{FinalNumUnmapped}
  >>>

  runtime{
    memory: "2 GB"
    disks: "local-disk 5 HDD"
  }
  output{
    Int NumberUnmappedReads = read_int("~{FinalNumUnmapped}")
    Int NumberMappedReads = read_int("~{FinalNumMapped}")
  }
}

