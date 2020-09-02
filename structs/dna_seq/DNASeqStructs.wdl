version 1.0

struct SampleAndUnmappedBams {
  String base_file_name
  String? final_gvcf_base_name
  Array[File] flowcell_unmapped_bams
  String sample_name
  String unmapped_bam_suffix
}

struct ReferenceFasta {
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File ref_alt
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac
}

struct DNASeqSingleSampleReferences {
  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File calling_interval_list

  ReferenceFasta reference_fasta

  Array[File] known_indels_sites_vcfs
  Array[File] known_indels_sites_indices

  File dbsnp_vcf
  File dbsnp_vcf_index

  File evaluation_interval_list

  File haplotype_database_file
}

struct VariantCallingScatterSettings {
   Int haplotype_scatter_count
   Int break_bands_at_multiples_of
}

struct ExomeGermlineSingleSampleOligos {
  File target_interval_list
  File bait_interval_list
  String bait_set_name
}

struct CrossSpeciesContaminationReferences {
  File filter_bwa_image
  File kmer_file
  File meats_bwa_image
  File meats_fasta
  File meats_fasta_dict
  File meats_taxonomy_file
  File microbe_bwa_image
  File microbe_fasta
  File microbe_fasta_dict
  File microbe_taxonomy_file
  File normalization_file
  File metrics_script_file
  Float score_min_identity
  Int reads_after_downsampling
}

struct PapiSettings {
  Int preemptible_tries
  Int agg_preemptible_tries
}
