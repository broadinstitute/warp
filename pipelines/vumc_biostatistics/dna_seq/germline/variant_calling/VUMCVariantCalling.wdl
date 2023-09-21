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

import "./VUMCHaplotypecallerReblockMoveResult.wdl" as Utils
import "../../../../../structs/dna_seq/DNASeqStructs.wdl"
import "../../../../../tasks/broad/BamProcessing.wdl" as Processing
import "../../../../broad/dna_seq/germline/variant_calling/VariantCalling.wdl" as ToGvcf

# WORKFLOW DEFINITION 
workflow VUMCVariantCalling {
  input {
    File input_bam
    File input_bam_index

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    String base_file_name
    String? final_gvcf_base_name

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    String gatk_path = "/gatk/gatk"
  
    Float? contamination

    Boolean use_spanning_event_genotyping = true

    String? project_id
    String? target_bucket
    String? genoset
    String? GRID
  }  
  
  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      use_spanning_event_genotyping = use_spanning_event_genotyping,
      calling_interval_list = references.calling_interval_list,
      evaluation_interval_list = references.evaluation_interval_list,
      haplotype_scatter_count = scatter_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = scatter_settings.break_bands_at_multiples_of,
      contamination = contamination,
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      ref_str = references.reference_fasta.ref_str,
      dbsnp_vcf = references.dbsnp_vcf,
      dbsnp_vcf_index = references.dbsnp_vcf_index,
      base_file_name = base_file_name,
      final_vcf_base_name = base_file_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries,
      run_dragen_mode_variant_calling = false,
      use_gatk3_haplotype_caller = false,
      use_dragen_hard_filtering = false
  }

  if(defined(target_bucket)){
    call Utils.MoveVcf {
      input:
        output_vcf = BamToGvcf.output_vcf,
        output_vcf_index = BamToGvcf.output_vcf_index,
        project_id = project_id,
        target_bucket = select_first([target_bucket]),
        genoset = select_first([genoset]),
        GRID = select_first([GRID])
    }
  }

  output {
    File vcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    File vcf_detail_metrics = BamToGvcf.vcf_detail_metrics
    File output_vcf = select_first([MoveVcf.target_output_vcf, BamToGvcf.output_vcf])
    File output_vcf_index = select_first([MoveVcf.target_output_vcf_index, BamToGvcf.output_vcf_index])
  }
}

