version 1.0

#WORKFLOW DEFINITION   
workflow VUMCCramQC2 {
  input {
    Array[File] input_crams
    String sample_name
    String samtools_docker = "https://hub.docker.com/u/mgibio/samtools-cwl"
    String samtools_path = "/samtools-cwl/"
    File reference_file
  }

  output {
    Int unmapped_reads = CRAM_unmapped_reads
    Int mapped_reads = CRAM_mapped_reads
  }

  call ValidateCRAM {
    input:
      input_crams = input_crams,
      sample_name = sample_name,
      samtools_docker = samtools_docker,
      samtools_path = samtools_path,
      reference_file = reference_file,
      
  }
}

# TASK DEFINITIONS
# Validate a cram using Picard ValidateSamFile
task ValidateCRAM {
  input {
    # Command parameters
    Array[File] input_crams
    String sample_name
    String samtools_path
    File reference_file

  
    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }
    
  Int disk_size = ceil(size(input_crams, "GB")) + addtional_disk_space_gb
  String NumUnmapped = "${sample_name}_Unmapped.txt"
  String NumMapped = "${sample_name}_Mapped.txt"

  command <<<
    echo "0" > ~{NumUnmapped}
    echo "0" > ~{NumMapped}

    for input_cram in ~{sep=" " input_crams}
    do
    
    ~{samtools_path} \
        samtools flagstat $input_cram |cut -f1 -d' '|head -n3|tail -n1 > ~{NumMapped}

    ~{samtools_path} \
        samtools view -c -T $reference_file $input_cram > ~{NumUnmapped}

    done
  >>>

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int NumberUnmappedReads = read_int("~{NumMapped}")
    Int NumberMappedReads = read_int("~{NumUnmapped}")
  }
}