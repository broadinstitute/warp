version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines tasks to use Dragen's DRAGstr approach to STR sequencing artifacts 
## Indel genotype priors in the DRAGEN-Gatk pipeline. 
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.


task ComposeSTRTableFile {
    
  input {
    File ref_fasta
    File ref_fasta_idx
    File ref_dict
    String docker = "us.gcr.io/broad-gatk/gatk:4.2.2.0"
    Int preemptible_tries = 3
  }

  String base_name = sub(basename(ref_fasta), ".fa(sta)$","")
  String str_table_file_name = base_name + ".str"
  Int disk_size_gb = ceil(size([ref_fasta, ref_fasta_idx], "GiB")) + 
                        21 # 1gb for the output file (usually far less than that) and 20 for the rest of the fs.

  command <<<

    gatk --java-options "-Xmx2g" \
      ComposeSTRTableFile \
        -R ~{ref_fasta} \
        -O ~{str_table_file_name}

  >>>

  runtime {
    docker: docker 
    disks: "local-disk " + disk_size_gb + " HDD"
    memory: "3 GiB"
    preemptible: preemptible_tries
  }

  output {
    File str_table_file = "~{str_table_file_name}"
  }
}

task CalibrateDragstrModel {

  input {
    File ref_fasta
    File ref_fasta_idx
    File ref_dict
    File str_table_file
    File alignment ## can handle cram or bam.
    File alignment_index
    String docker = "us.gcr.io/broad-gatk/gatk:4.2.2.0"
    Int preemptible_tries = 3
    Int threads = 4
    Int? mem_gb
    Boolean use_ssd = true
  }

  # If CRAM, restrict threads to a maximum of 4
  Boolean is_cram = sub(alignment, "\\.cram$", "") != "" + alignment
  Int java_threads = if (threads < 1 ) then 1 
                else if (is_cram && threads > 4) then 4 # more than 4 threads in cram is probrably contra-productive.
                else threads

  String base_name = basename(alignment)
  String out_file_name = base_name + ".dragstr"
  Int disk_size_gb = ceil(size([ref_fasta, ref_fasta_idx, alignment, str_table_file], "GiB")) + 
                        20 # 20 for the rest of the fs.

  String parallel_args  = if (java_threads <= 1) then "" else "--threads " + threads
  
  # If the input is a CRAM we need an additional 500MB of memory per thread
  Int recommended_memory_gb = ceil(2 + (if (is_cram) then 0.5 else 0.1) * java_threads)
  Int runtime_memory_gb = select_first([mem_gb, recommended_memory_gb])
  Int java_memory_gb = if (runtime_memory_gb < 2) then 1 else runtime_memory_gb - 1

  command <<<
    set -x
    gatk --java-options "-Xmx~{java_memory_gb}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Dsamjdk.reference_fasta=~{ref_fasta}" \
      CalibrateDragstrModel \
        -R ~{ref_fasta} \
        -I ~{alignment} \
        -str ~{str_table_file} \
        -O ~{out_file_name} \
        ~{parallel_args}

  >>>

  runtime {
     docker: docker
     disks: "local-disk " + disk_size_gb + (if use_ssd then " SSD" else " HDD")
     memory: (java_memory_gb + 1) + " GiB"
     preemptible: preemptible_tries
     cpu: java_threads
  }

  output {
    File dragstr_model = "~{out_file_name}"
  }
}
