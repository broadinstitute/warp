version 1.0

## Copyright Broad Institute/VUMC, 2018/2022
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in FASTQ format
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
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

import "../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow VUMCFastqToUnmappedCram {

  String pipeline_version = "3.1.10"

  input {
    # Optional for VUMC pipeline
    String sample_name 
    String fastq_1 
    String fastq_2 
    String readgroup_name 
    String? library_name 
    String? platform_unit 
    String? run_date 
    String? platform_name 
    String? sequencing_center 
  }

  # Convert pair of FASTQs to uCRAM
  call PairedFastQsToUnmappedCram {
    input:
      sample_name = sample_name,
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,
      readgroup_name = readgroup_name,
      library_name = library_name,
      platform_unit = platform_unit,
      run_date = run_date,
      platform_name = platform_name,
      sequencing_center = sequencing_center,
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_unmapped_cram = PairedFastQsToUnmappedCram.output_unmapped_cram
  }
  meta {
    allowNestedInputs: true
  }
}

# Convert a pair of FASTQs to uCRAM
task PairedFastQsToUnmappedCram {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String? library_name 
    String? platform_unit 
    String? run_date 
    String? platform_name 
    String? sequencing_center 

    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 7
    Int preemptible_attempts = 3

    # The cram file should be smaller than the original FASTQ files.
    # In these cases we need to account for the input (1) and the output(1).
    Float disk_multiplier = 1

    String docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
  }
  Int command_mem_gb = machine_mem_gb - 1
  Float fastq_size = size(fastq_1, "GB") + size(fastq_2, "GB")
  Int disk_space_gb = ceil(fastq_size + (fastq_size * disk_multiplier ) + addtional_disk_space_gb)
  command {
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
      FastqToSam \
      --FASTQ ~{fastq_1} \
      --FASTQ2 ~{fastq_2} \
      --OUTPUT ~{sample_name}.unmapped.cram \
      --SAMPLE_NAME ~{sample_name} \
      ~{"--LIBRARY_NAME " + library_name} \
      ~{"--PLATFORM_UNIT " + platform_unit} \
      ~{"--RUN_DATE " + run_date} \
      ~{"--PLATFORM " + platform_name} \
      ~{"--SEQUENCING_CENTER " + sequencing_center} \
      --READ_GROUP_NAME ~{readgroup_name}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
  }
  output {
    File output_unmapped_cram = "~{sample_name}.unmapped.cram"
  }
}
