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

struct RuntimeOptions {
  Int? preemptible_tries
  Int? evaluation_preemptible_tries
  String monitoring_script
  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size
  Int? additional_metrics_disk # will be added to increase_disk_size
  Boolean? no_address_override
  # When running on Terra, use workspace.name as this input to ensure that all tasks will only cache hit to runs in your
  # own workspace. This will prevent call caching from failing with "Cache Miss (10 failed copy attempts)". Outside of
  # Terra this can be left as the default empty String. This dummy input is only needed for tasks that have no inputs
  # specific to the sample being run (such as GetBwaVersion which does not take in any sample data).
  String dummy_input_for_call_caching
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

struct EnvironmentVersions {
  String broad_gatk_docker
  String crammer_docker
  String jb_gatk_docker
  String gatk_markduplicates_docker
  String jukebox_vc_docker
  String gitc_docker
  String? gitc_path_override
  File? picard_jar_override
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