version 1.0

#WORKFLOW DEFINITION
workflow VUMCCramQC2 {
  input {
    Array[File] input_crams
    String sample_name
    String samtools_docker = "staphb/samtools:latest"
    File reference_file
  }

  
 scatter (input_cram in input_crams){
  call CountCRAM {
    input:
      input_cram = input_cram,
      docker = samtools_docker,
      reference_file = reference_file
  }
 }

 call SumUp {
  input: 
   sample_name = sample_name,
   mapped_reads = CountCRAM.mapped_reads,
   unmapped_reads = CountCRAM.unmapped_reads,
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
    File reference_file

    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }

  Int disk_size = ceil(size(input_cram, "GB")) + addtional_disk_space_gb
  String NumUnmapped = "Unmapped.txt"
  String NumMapped = "Mapped.txt"

  command <<<
    samtools flagstat $input_cram |grep "mapped (" |cut -f1 -d' ' >> ~{NumMapped}
    samtools view -c -f4 -T ~{reference_file} $input_cram >> ~{NumUnmapped}
  >>>

  runtime{
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output{
   File unmapped_reads = "~{NumUnmapped}"
   File mapped_reads = "~{NumMapped}"
  }
  
}

task SumUp{
  input{
    # Command parameters
     String sample_name
     Array[File] unmapped_reads
     Array[File] mapped_reads

    # Runtime parameters
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }

  String FinalNumUnmapped = "${sample_name}_final_Unmapped.txt"
  String FinalNumMapped = "${sample_name}_final_Mapped.txt"

  command <<<
    cat ~{mapped_reads} | awk '{ sum += $1 } END { print sum }' > ~{FinalNumMapped}
    cat ~{unmapped_reads} | awk '{ sum += $1 } END { print sum }' > ~{FinalNumUnmapped}
  >>>

  runtime{
    memory: machine_mem_gb + " GB"
    disks: "local-disk "
  }
  output{
    Int NumberUnmappedReads = read_int("~{FinalNumUnmapped}")
    Int NumberMappedReads = read_int("~{FinalNumMapped}")
  }
}

