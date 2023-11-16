version 1.0

task RevertSamSingle {
  input {
    #Command parameters
    File input_cram
    String input_cram_suffix

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    #Runtime parameters
    Int machine_mem_gb = 10
    Int preemptible_attempts = 3

    Float disk_multiplier = 3
    Int additional_disk_size = 10

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  String output_file = basename(input_cram, input_cram_suffix) + ".unsorted." + input_cram_suffix

  Int disk_size = ceil(size(input_cram, "GB") * disk_multiplier + additional_disk_size)
  Int command_mem_gb = machine_mem_gb - 1    ####Needs to occur after machine_mem_gb is set 

  command <<<
 
    gatk --java-options "-Xmx~{command_mem_gb}g" \
    RevertSam \
    --INPUT ~{input_cram} \
    --REFERENCE_SEQUENCE ~{ref_fasta} \
    --OUTPUT ~{output_file} \
    --OUTPUT_BY_READGROUP false \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER coordinate
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " HDD"
    memory: machine_mem_gb + " GB"
    preemptible: preemptible_attempts
  }
  output {
    File unmapped_crams = "~{output_file}"
  }
}

task RevertSamByReadGroup {
  input {
    #Command parameters
    File input_cram
    String input_cram_suffix

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Boolean restore_hardclips = true

    #Runtime parameters
    Int machine_mem_gb = 10
    Int preemptible_attempts = 3

    Float disk_multiplier = 6
    Int additional_disk_size = 20

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  Int disk_size = ceil(size(input_cram, "GB") * disk_multiplier + additional_disk_size)
  Int command_mem_gb = machine_mem_gb - 1    ####Needs to occur after machine_mem_gb is set 

  command <<<
 
    gatk --java-options "-Xmx~{command_mem_gb}g" \
    RevertSam \
    --INPUT ~{input_cram} \
    --REFERENCE_SEQUENCE ~{ref_fasta} \
    --OUTPUT ./ \
    --RESTORE_HARDCLIPS ~{restore_hardclips} \
    --OUTPUT_BY_READGROUP true \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER coordinate
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " HDD"
    memory: machine_mem_gb + " GB"
    preemptible: preemptible_attempts
  }
  output {
    Array[File] unmapped_crams = glob("*" + input_cram_suffix)
  }
}

task SortSam {
  input {
    #Command parameters
    File input_cram
    String sorted_bam_name

    #Runtime parameters
    Int machine_mem_gb = 10
    Int preemptible_attempts = 3

    Float disk_multiplier = 6
    Int additional_disk_size = 10
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  Int disk_size = ceil(size(input_cram, "GB") * disk_multiplier + additional_disk_size)
  Int command_mem_gb = machine_mem_gb - 1    ####Needs to occur after machine_mem_gb is set 

  command <<<
    gatk --java-options "-Xmx~{command_mem_gb}g" \
    SortSam \
    --INPUT ~{input_cram} \
    --OUTPUT ~{sorted_bam_name} \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 1000000
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " HDD"
    memory: machine_mem_gb + " GB"
    preemptible: preemptible_attempts
  }
  output {
    File sorted_bam = "~{sorted_bam_name}"
  }
}

