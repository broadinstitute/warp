version 1.0

task VerifyPipelineInputs {
  input {
    File? bam
    File? r1_fastq
    File? r2_fastq
  
    # Only needed if inputs are fastq rather than ubam
	  String? library_name
	  String? platform
	  String? platform_unit
	  String? read_group_name
	  String? sequencing_center = "BI"

    String docker = "python:3.7.2"
    Int cpu = 1
    Int memory_mb = 1000
    Int disk_size_gb = 1
  }

  command <<<
  set -e
  python3 <<CODE

  fastq_flag = 0
  bam = "~{bam}"
  r1_fastq = "~{r1_fastq}"
  r2_fastq = "~{r2_fastq}"
  library_name = "~{library_name}"
	platform = "~{platform}"
	platform_unit = "~{platform_unit}"
	read_group_name = "~{read_group_name}"
	sequencing_center = "~{sequencing_center}"

  if bam and fastq1 or fastq2:
    raise ValueError("Bam and Fastq files cannot both be defined as input")
  
  if not bam and not r1_fastq and not r2_fastq:
    raise ValueError("Either Bam or Fastq must be defined")

  if not bam:
    if r1_fastq and r2_fastq:
      if not library_name or not platform or not platform_unit or not read_group_name or not sequencing_center:
        raise ValueError("r1_fastq and r2_fastq defined, must have: library_name, platform, platform_unit, read_group_name, sequencing_center")
      else:
        fastq_flag += 1
    else
      raise ValueError("r1_fastq and r2_fastq must both be defined")

  
  with open("output.txt", as "w") as f:
    if fastq_flag == 1:
      f.write("true")
    else: # Remaining case is that only bam is defined
      f.write("false")

  CODE
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    Boolean fastq_run = read_boolean("output.txt")
  }

}