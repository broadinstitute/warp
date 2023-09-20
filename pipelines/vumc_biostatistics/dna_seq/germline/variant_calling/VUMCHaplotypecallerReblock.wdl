version 1.0

## Copyright Broad Institute, 2019
## 
## The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool
## from GATK4 in GVCF mode on a single sample according to GATK Best Practices.
## When executed the workflow scatters the HaplotypeCaller tool over a sample
## using an intervals list file. The output file produced will be a
## single gvcf file which can be used by the joint-discovery workflow.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##
## Cromwell version support 
## - Successfully tested on v53
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

import "../../../../../tasks/broad/Utilities.wdl" as Utilities
import "../../../../broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl" as Reblock
import "./VUMCHaplotypecallerReblockMoveResult.wdl" as Utils

# WORKFLOW DEFINITION 
workflow VUMCHaplotypecallerReblock {
  input {
    File input_bam
    File input_bam_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File scattered_calling_intervals_list

    String? project_id
    String? target_bucket
    String? genoset
    String? GRID

    Boolean make_gvcf = true
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    String gatk_path = "/gatk/gatk"
  }  

  # if(defined(target_bucket)){
  #   if(!defined(genoset)){
  #     call Utilities.ErrorWithMessage as NoGenosetError {
  #       input:
  #         message = "genoset is missing when target bucket is set."
  #     }
  #   }
  #   if(!defined(GRID)){
  #     call Utilities.ErrorWithMessage as GRIDError {
  #       input:
  #         message = "GRID is missing when target bucket is set."
  #     }
  #   }
  # }

  Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

  #is the input a cram file?
  Boolean is_cram = sub(basename(input_bam), ".*\\.", "") == "cram"

  String sample_basename = if is_cram then  basename(input_bam, ".cram") else basename(input_bam, ".bam")
  String vcf_basename = sample_basename
  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_filename = vcf_basename + output_suffix

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = length(scattered_calling_intervals) - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over grouped calling intervals
  scatter (interval_file in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        interval_list = interval_file,
        output_filename = output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        make_gvcf = make_gvcf,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      output_filename = output_filename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  call Reblock.ReblockGVCF as Reblock {
    input:
      gvcf = MergeGVCFs.output_vcf,
      gvcf_index = MergeGVCFs.output_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      calling_interval_list = scattered_calling_intervals_list
  }

  if(defined(target_bucket)){
    call Utils.MoveVcf {
      input:
        output_vcf = Reblock.output_vcf,
        output_vcf_index = Reblock.output_vcf_index,
        project_id = project_id,
        target_bucket = select_first([target_bucket]),
        genoset = select_first([genoset]),
        GRID = select_first([GRID])
    }
  }

  output {
    File output_vcf = select_first([MoveVcf.target_output_vcf, Reblock.output_vcf])
    File output_vcf_index = select_first([MoveVcf.target_output_vcf_index, Reblock.output_vcf_index])
  }
}

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    File interval_list
    String output_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Int hc_scatter

    String? gcs_project_for_requester_pays

    String gatk_path
    String? java_options

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
  }

  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"]) 

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

  String vcf_basename = if make_gvcf then  basename(output_filename, ".gvcf") else basename(output_filename, ".vcf")

  parameter_meta {
    input_bam: {
      description: "a bam/cram file",
      localization_optional: true
    }
    input_bam_index: {
      description: "an index file for the bam/cram input",
      localization_optional: true
    }
  }
  command {
    set -e
  
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_filename} \
      -contamination ~{default="0" contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""} 
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename

    String gatk_path

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1
  
  command {
  set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}

