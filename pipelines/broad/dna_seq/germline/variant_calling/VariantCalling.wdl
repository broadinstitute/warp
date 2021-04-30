version 1.0

import "../../../../../tasks/broad/GermlineVariantDiscovery.wdl" as Calling
import "../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../tasks/broad/Utilities.wdl" as Utils
import "../../../../../tasks/broad/BamProcessing.wdl" as BamProcessing

workflow VariantCalling {

  String pipeline_version = "1.0.0"

  input {
    File calling_interval_list
    File evaluation_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Float? contamination
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp_vcf
    File dbsnp_vcf_index
    String base_file_name
    String final_vcf_base_name
    Int agg_preemptible_tries
    Boolean make_gvcf = true
    Boolean make_bamout = false
    Boolean use_gatk3_haplotype_caller = false
  }

  parameter_meta {
    make_bamout: "For CNNScoreVariants to run with a 2D model, a bamout must be created by HaplotypeCaller. The bamout is a bam containing information on how HaplotypeCaller remapped reads while it was calling variants. See https://gatkforums.broadinstitute.org/gatk/discussion/5484/howto-generate-a-bamout-file-showing-how-haplotypecaller-has-remapped-sequence-reads for more details."
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utils.ScatterIntervalList as ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
  scatter (scattered_interval_list in ScatterIntervalList.out) {

    if (use_gatk3_haplotype_caller) {
      call Calling.HaplotypeCaller_GATK35_GVCF as HaplotypeCallerGATK3 {
        input:
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = scattered_interval_list,
          gvcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          contamination = contamination,
          preemptible_tries = agg_preemptible_tries,
          hc_scatter = hc_divisor
      }
    }

    if (!use_gatk3_haplotype_caller) {

      # Generate GVCF by interval
      call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
        input:
          contamination = contamination,
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = scattered_interval_list,
          vcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          hc_scatter = hc_divisor,
          make_gvcf = make_gvcf,
          make_bamout = make_bamout,
          preemptible_tries = agg_preemptible_tries
       }

      # If bamout files were created, we need to sort and gather them into one bamout
      if (make_bamout) {
        call BamProcessing.SortSam as SortBamout {
          input:
            input_bam = HaplotypeCallerGATK4.bamout,
            output_bam_basename = final_vcf_base_name,
            preemptible_tries = agg_preemptible_tries,
            compression_level = 2
        }
      }
    }

    File vcfs_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf, HaplotypeCallerGATK4.output_vcf])
    File vcf_indices_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf_index, HaplotypeCallerGATK4.output_vcf_index])
  }

  # Combine by-interval (g)VCFs into a single sample (g)VCF file
  String merge_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  call Calling.MergeVCFs as MergeVCFs {
    input:
      input_vcfs = vcfs_to_merge,
      input_vcfs_indexes = vcf_indices_to_merge,
      output_vcf_name = final_vcf_base_name + merge_suffix,
      preemptible_tries = agg_preemptible_tries
  }

  if (make_bamout) {
    call MergeBamouts {
      input:
        bams = select_all(SortBamout.output_bam),
        output_base_name = final_vcf_base_name
    }
  }

  # Validate the (g)VCF output of HaplotypeCaller
  call QC.ValidateVCF as ValidateVCF {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      calling_interval_list = calling_interval_list,
      is_gvcf = make_gvcf,
      preemptible_tries = agg_preemptible_tries
  }

  # QC the (g)VCF
  call QC.CollectVariantCallingMetrics as CollectVariantCallingMetrics {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      metrics_basename = final_vcf_base_name,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      ref_dict = ref_dict,
      evaluation_interval_list = evaluation_interval_list,
      is_gvcf = make_gvcf,
      preemptible_tries = agg_preemptible_tries
  }

  output {
    File vcf_summary_metrics = CollectVariantCallingMetrics.summary_metrics
    File vcf_detail_metrics = CollectVariantCallingMetrics.detail_metrics
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
    File? bamout = MergeBamouts.output_bam
    File? bamout_index = MergeBamouts.output_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}

# This task is here because merging bamout files using Picard produces an error.
task MergeBamouts {

  input {
    Array[File] bams
    String output_base_name
  }

  Int disk_size = ceil(size(bams, "GiB") * 2) + 10

  command {
    samtools merge ~{output_base_name}.bam ~{sep=" " bams}
    samtools index ~{output_base_name}.bam
    mv ~{output_base_name}.bam.bai ~{output_base_name}.bai
  }

  output {
    File output_bam = "~{output_base_name}.bam"
    File output_bam_index = "~{output_base_name}.bai"
  }

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    memory: "4 GiB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 3
    cpu: 1
  }
}