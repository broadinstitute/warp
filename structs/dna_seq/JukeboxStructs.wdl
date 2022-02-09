version 1.0

struct SampleInputs {
  Array[File] input_cram_bam_list
  String base_file_name
  String? override_input_ending # For drs where there is no extension. Should be "is_cram" or "is_bam"
}

struct ContaminationSites {
  String contamination_sites_path
  File contamination_sites_vcf
  File contamination_sites_vcf_index
}

struct References {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
}

struct AlignmentReferences {
  References references
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
}

struct VariantCallingSettings {
  File wgs_calling_interval_list
  Int break_bands_at_multiples_of
  Int haplotype_scatter_count
}

struct ExtraArgs {
  Int? reads_per_split
  String? hc_extra_args
  String? mark_duplicates_extra_args
  Float rsq_threshold
}

struct VcfPostProcessing {
  Array[File] annotation_intervals
  File? filtering_model_no_gt
  File af_only_gnomad
  File af_only_gnomad_index
  Boolean filter_cg_insertions
  File? filtering_blacklist_file
  File? training_blacklist_file
  Int? exome_weight
  String? exome_weight_annotation
  File? interval_list_override
  File runs_file
  String? filtering_model_no_gt_name_override
  String? filtering_model_with_gt_name_override
  Float max_duplication_in_reasonable_sample
  Float max_chimerism_in_reasonable_sample
  Boolean? make_gvcf_override
  Boolean? merge_bam_file_override
  File ref_dbsnp
  File ref_dbsnp_index
  File wgs_coverage_interval_list
}