version 1.0

import "../../../somatic/single_sample/ugwgs/UltimaGenomicsWholeGenomeCramOnly.wdl" as UltimaGenomicsWholeGenomeCramOnly
import "../../../../../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineTasks.wdl" as Tasks
import "../../../../../../tasks/wdl/Utilities.wdl" as Utilities
import "../../../../../../tasks/wdl/GermlineVariantDiscovery.wdl" as VariantDiscoverTasks
import "../../../../../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates.wdl" as UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates
import "../../../../../../tasks/wdl/InternalTasks.wdl" as InternalTasks
import "../../../../../../tasks/wdl/Qc.wdl" as QC
import "../../../../../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineQC.wdl" as UltimaGenomicsWholeGenomeGermlineQC
import "../../../../../../structs/dna_seq/UltimaGenomicsWholeGenomeGermlineStructs.wdl" as Structs
import "../../joint_genotyping/reblocking/ReblockGVCF.wdl" as ReblockGVCF


workflow UltimaGenomicsWholeGenomeGermline {
  input {

    ContaminationSites contamination_sites
    AlignmentReferences alignment_references
    VariantCallingSettings variant_calling_settings
    VcfPostProcessing vcf_post_processing

    # Sample Information
    Array[File]? input_cram_list
    Array[File]? input_bam_list
    String base_file_name

    Float rsq_threshold = 1.0
    Boolean merge_bam_file = true
    Boolean make_haplotype_bam = false
    Int reads_per_split = 20000000
    String filtering_model_no_gt_name = "rf_model_ignore_gt_incl_hpol_runs"
  }

  meta {
    allowNestedInputs: true
  }

  parameter_meta {
    contamination_sites: "Struct containing files for contamination estimation"
    alignment_references: "Struct containing reference files for alignment with BWA mem"
    variant_calling_settings: "Struct containing reference files for variant calling with HaplotypeCaller"
    vcf_post_processing: "Struct containing reference files for VCF post-processing: annotation and filtering"
    input_cram_list: "Array of CRAM files to be used as workflow input. Must be specified if `input_bam_list` is not provided"
    input_bam_list: "Array of unmapped BAM files to be used as workflow input. Must be specified if `input_cram_list` is not provided"
    base_file_name: "Base name for each of the output files."
    rsq_threshold: "Threshold for a read quality metric that is produced by the sequencing platform"
    merge_bam_file: "Boolean indicating if by-interval bamout files from HaplotypeCaller should be merged into a single BAM"
    reads_per_split: "Number of reads by which to split the CRAM prior to alignment"
    filtering_model_no_gt_name: "String describing the optional filtering model; default set to rf_model_ignore_gt_incl_hpol_runs"
  }

  String pipeline_version = "1.2.1"


  References references = alignment_references.references

  call UltimaGenomicsWholeGenomeCramOnly.UltimaGenomicsWholeGenomeCramOnly {
    input:
      contamination_sites = contamination_sites,
      alignment_references = alignment_references,
      input_cram_list = input_cram_list,
      input_bam_list = input_bam_list,
      base_file_name = base_file_name,
      rsq_threshold = rsq_threshold,
      reads_per_split = reads_per_split,
      vcf_post_processing = vcf_post_processing,
      save_bam_file = true
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utilities.ScatterIntervalList {
    input:
      interval_list               = variant_calling_settings.wgs_calling_interval_list,
      scatter_count               = variant_calling_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = variant_calling_settings.break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by factor 2 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int hc_divisor = ScatterIntervalList.interval_count / 2

  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate VCF by interval
    call Tasks.HaplotypeCaller as HaplotypeCaller {
      input:
        input_bam_list  = [select_first([UltimaGenomicsWholeGenomeCramOnly.output_bam])],
        input_bam_index_list = [select_first([UltimaGenomicsWholeGenomeCramOnly.output_bam_index])],
        interval_list   = ScatterIntervalList.out[index],
        vcf_basename    = UltimaGenomicsWholeGenomeCramOnly.output_safe_name,
        references      = references,
        make_bamout     = make_haplotype_bam
    }
  }

  # Combine by-interval VCFs into a single sample  file
  call VariantDiscoverTasks.MergeVCFs {
    input:
      input_vcfs          = HaplotypeCaller.output_vcf,
      input_vcfs_indexes  = HaplotypeCaller.output_vcf_index,
      output_vcf_name     = UltimaGenomicsWholeGenomeCramOnly.output_safe_name + ".g.vcf.gz",
  }

  # Combine by-interval BAMs into a single sample file
  if(merge_bam_file && make_haplotype_bam) {
      call Tasks.MergeBams {
        input:
          input_bams      = HaplotypeCaller.bamout,
          output_bam_name = UltimaGenomicsWholeGenomeCramOnly.output_safe_name + ".bam",
      }
  }

  call Tasks.ConvertGVCFtoVCF {
    input:
      input_gvcf         = MergeVCFs.output_vcf,
      input_gvcf_index   = MergeVCFs.output_vcf_index,
      output_vcf_name    = UltimaGenomicsWholeGenomeCramOnly.output_safe_name + '.vcf.gz',
      references         = references
  }

  # VCF post-processings
  call Tasks.AnnotateVCF {
    input :
      input_vcf               = ConvertGVCFtoVCF.output_vcf,
      input_vcf_index         = ConvertGVCFtoVCF.output_vcf_index,
      references              = references,
      reference_dbsnp         = vcf_post_processing.ref_dbsnp,
      reference_dbsnp_index   = vcf_post_processing.ref_dbsnp_index,
      flow_order              = UltimaGenomicsWholeGenomeCramOnly.flow_order,
      final_vcf_base_name     = UltimaGenomicsWholeGenomeCramOnly.output_safe_name
  }

  call Tasks.AddIntervalAnnotationsToVCF {
    input:
      input_vcf             = AnnotateVCF.output_vcf_annotated,
      input_vcf_index       = AnnotateVCF.output_vcf_annotated_index,
      final_vcf_base_name   = UltimaGenomicsWholeGenomeCramOnly.output_safe_name,
      annotation_intervals  = vcf_post_processing.annotation_intervals
  }

  call Tasks.TrainModel {
    input:
      input_file                = AddIntervalAnnotationsToVCF.output_vcf,
      input_file_index          = AddIntervalAnnotationsToVCF.output_vcf_index,
      input_vcf_name            = UltimaGenomicsWholeGenomeCramOnly.output_safe_name,
      blocklist_file            = vcf_post_processing.training_blocklist_file,
      ref_fasta                 = references.ref_fasta,
      ref_index                 = references.ref_fasta_index,
      runs_file                 = vcf_post_processing.runs_file,
      apply_model               = filtering_model_no_gt_name,
      annotation_intervals      = vcf_post_processing.annotation_intervals,
      exome_weight              = vcf_post_processing.exome_weight,
      exome_weight_annotation   = vcf_post_processing.exome_weight_annotation
  }


  call Tasks.AnnotateVCF_AF {
    input :
      input_vcf             = AddIntervalAnnotationsToVCF.output_vcf,
      input_vcf_index       = AddIntervalAnnotationsToVCF.output_vcf_index,
      af_only_gnomad        = vcf_post_processing.af_only_gnomad,
      af_only_gnomad_index  = vcf_post_processing.af_only_gnomad_index,
      final_vcf_base_name   = UltimaGenomicsWholeGenomeCramOnly.output_safe_name
  }

  call Tasks.FilterVCF {
    input:
      input_vcf               = AnnotateVCF_AF.output_vcf_annotated,
      input_model             = select_first([vcf_post_processing.filtering_model_no_gt,TrainModel.model_pkl]),
      runs_file               = vcf_post_processing.runs_file,
      references              = references,
      model_name              = filtering_model_no_gt_name,
      filter_cg_insertions    = vcf_post_processing.filter_cg_insertions,
      final_vcf_base_name     = UltimaGenomicsWholeGenomeCramOnly.output_safe_name,
      flow_order              = UltimaGenomicsWholeGenomeCramOnly.flow_order,
      annotation_intervals    = vcf_post_processing.annotation_intervals
  }

  call Tasks.MoveAnnotationsToGvcf {
    input:
      filtered_vcf        = FilterVCF.output_vcf_filtered,
      filtered_vcf_index  = FilterVCF.output_vcf_filtered_index,
      gvcf                = MergeVCFs.output_vcf,
      gvcf_index          = MergeVCFs.output_vcf_index
  }

  call ReblockGVCF.ReblockGVCF {
    input:
      gvcf = MoveAnnotationsToGvcf.output_gvcf,
      gvcf_index = MoveAnnotationsToGvcf.output_gvcf_index,
      calling_interval_list = variant_calling_settings.wgs_calling_interval_list,
      ref_dict = alignment_references.references.ref_dict,
      ref_fasta = alignment_references.references.ref_fasta,
      ref_fasta_index = alignment_references.references.ref_fasta_index,
      tree_score_cutoff = vcf_post_processing.remove_low_tree_score_sites_cutoff,
      annotations_to_keep_command = vcf_post_processing.annotations_to_keep_command_for_reblocking,
      cloud_provider = "gcp"
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_gvcf = ReblockGVCF.reblocked_gvcf
    File output_gvcf_index = ReblockGVCF.reblocked_gvcf_index
    File output_vcf = ConvertGVCFtoVCF.output_vcf
    File output_vcf_index = ConvertGVCFtoVCF.output_vcf_index

    #MERGE bam file
    File? haplotype_bam = MergeBams.output_bam
    File? haplotype_bam_index = MergeBams.output_bam_index

    File output_cram = UltimaGenomicsWholeGenomeCramOnly.output_cram
    File output_cram_index = UltimaGenomicsWholeGenomeCramOnly.output_cram_index
    File output_cram_md5 = UltimaGenomicsWholeGenomeCramOnly.output_cram_md5

    File selfSM = UltimaGenomicsWholeGenomeCramOnly.selfSM
    Float contamination = UltimaGenomicsWholeGenomeCramOnly.contamination

    # VCF post-processing
    File filtered_vcf = FilterVCF.output_vcf_filtered
    File filtered_vcf_index = FilterVCF.output_vcf_filtered_index

    # STATISTIC COLLECTION
    File quality_yield_metrics = UltimaGenomicsWholeGenomeCramOnly.quality_yield_metrics
    File wgs_metrics = UltimaGenomicsWholeGenomeCramOnly.wgs_metrics
    File raw_wgs_metrics = UltimaGenomicsWholeGenomeCramOnly.raw_wgs_metrics
    File duplicate_metrics = UltimaGenomicsWholeGenomeCramOnly.duplicate_metrics
    File agg_alignment_summary_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_alignment_summary_metrics
    File? agg_alignment_summary_pdf = UltimaGenomicsWholeGenomeCramOnly.agg_alignment_summary_pdf
    File agg_gc_bias_detail_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_gc_bias_summary_metrics
    File agg_quality_distribution_pdf = UltimaGenomicsWholeGenomeCramOnly.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = UltimaGenomicsWholeGenomeCramOnly.agg_quality_distribution_metrics
    Float duplication_rate = UltimaGenomicsWholeGenomeCramOnly.duplication_rate
    Float chimerism_rate = UltimaGenomicsWholeGenomeCramOnly.chimerism_rate
    Boolean is_outlier_data = UltimaGenomicsWholeGenomeCramOnly.is_outlier_data

    String sample_name = UltimaGenomicsWholeGenomeCramOnly.sample_name
    String flow_order = UltimaGenomicsWholeGenomeCramOnly.flow_order
    String barcode = UltimaGenomicsWholeGenomeCramOnly.barcode
    String id = UltimaGenomicsWholeGenomeCramOnly.id
  }

}
