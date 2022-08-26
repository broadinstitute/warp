version 1.0
##Copyright Broad Institute, 2018, VUMC/Biosttistics, 2022
## 
## This WDL converts multiple paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (multiple files per orientation)
## - The following metada descriptors per sample:
##  - readgroup
##  - sample_name
##  - library_name
##  - platform_unit
##  - run_date
##  - platform_name
##  - sequecing_center
##
## Outputs :
## - One unmapped bam
##
## Cromwell version support 
## - Successfully tested on v47
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow ConvertMultiplePairedFastQsToUnmappedBamWf {
  input {
    String sample_name 
    Array[File] fastq_1_list 
    Array[File] fastq_2_list 
    String readgroup_name 
    String library_name 
    String platform_unit 
    String run_date 
    String platform_name 
    String sequencing_center 

    Boolean make_fofn = false

    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
    
    # Sometimes the output is larger than the input, or a task can spill to disk.
    # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
    Float disk_multiplier = 2.5
  }

  String ubam_list_name = sample_name

  Array[Pair[File,File]] fastq_pairs = zip(fastq_1_list, fastq_2_list)

  scatter (fastq in fastq_pairs) {
    # Convert pair of FASTQs to uBAM
    call PairedFastQsToUnmappedBAM {
      input:
        sample_name = sample_name,
        fastq_1 = fastq.left,
        fastq_2 = fastq.right,
        readgroup_name = readgroup_name,
        library_name = library_name,
        platform_unit = platform_unit,
        run_date = run_date,
        platform_name = platform_name,
        sequencing_center = sequencing_center,
        gatk_path = gatk_path,
        docker = gatk_docker,
        disk_multiplier = disk_multiplier
    }

    Float current_mapped_size = size(PairedFastQsToUnmappedBAM.output_unmapped_bam, "GiB")    
  }

  call SumFloats as SumUnmappedBamSizes {
    input:
      sizes = current_mapped_size,
      preemptible_tries = preemptible_tries
  }

  call GatherUnmappedBamFiles {
    input:
      input_bams = PairedFastQsToUnmappedBAM.output_unmapped_bam,
      total_input_size = SumUnmappedBamSizes.total_size,
      output_bam_basename = sample_name,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }  

  # Outputs that will be retained when execution is complete
  output {
    File output_unmapped_bam = GatherUnmappedBamFiles.output_unmapped_bam
  }
}

# TASK DEFINITIONS

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center
    String gatk_path

    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 7
    Int preemptible_attempts = 3
    String docker
    Float disk_multiplier
  }
    Int command_mem_gb = machine_mem_gb - 1
    Float fastq_size = size(fastq_1, "GB") + size(fastq_2, "GB")
    Int disk_space_gb = ceil(fastq_size + (fastq_size * disk_multiplier ) + addtional_disk_space_gb)
  command {
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM_UNIT ~{platform_unit} \
    --RUN_DATE ~{run_date} \
    --PLATFORM ~{platform_name} \
    --SEQUENCING_CENTER ~{sequencing_center} 
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
  }
  output {
    File output_unmapped_bam = "~{readgroup_name}.unmapped.bam"
  }
}

# Calculates sum of a list of floats
task SumFloats {
  input {
    Array[Float] sizes
    Int preemptible_tries
  }

  command <<<
    python3 -c 'print(~{sep="+" sizes})'
  >>>
  output {
    Float total_size = read_float(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: preemptible_tries
  }
}

task GatherUnmappedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_unmapped_bam = "~{output_bam_basename}.unmapped.bam"
  }
}

