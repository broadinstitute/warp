version 1.0

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String? output_bam_basename

    String library_name 
    String platform_unit 
    String platform_name 
    String? run_date 
    String? sequencing_center 

    # Runtime parameters
    Int addtional_disk_space_gb = 100
    Int machine_mem_gb = 7
    Int preemptible_attempts = 3

    # Sometimes the output is larger than the input, or a task can spill to disk.
    # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
    Float disk_multiplier = 2.5

    String docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
  }

  String output_file = select_first([output_bam_basename, sample_name]) + ".unmapped.bam"
  Int command_mem_gb = machine_mem_gb - 1
  Float fastq_size = size(fastq_1, "GB") + size(fastq_2, "GB")
  Int disk_space_gb = ceil(fastq_size + (fastq_size * disk_multiplier ) + addtional_disk_space_gb)
  command <<<
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
      FastqToSam \
      --FASTQ ~{fastq_1} \
      --FASTQ2 ~{fastq_2} \
      --OUTPUT ~{output_file} \
      --SAMPLE_NAME ~{readgroup_name} \
      ~{"--LIBRARY_NAME " + library_name} \
      ~{"--PLATFORM_UNIT " + platform_unit} \
      ~{"--RUN_DATE " + run_date} \
      ~{"--PLATFORM " + platform_name} \
      ~{"--SEQUENCING_CENTER " + sequencing_center} \
      --READ_GROUP_NAME ~{readgroup_name}
  >>>
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
  }
  output {
    File output_unmapped_bam = "~{output_file}"
  }
}
