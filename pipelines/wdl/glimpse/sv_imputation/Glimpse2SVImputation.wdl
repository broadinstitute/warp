version 1.0

import "./PreprocessPLsGVCF.wdl" as PreprocessPLsGVCF
import "./Glimpse2SVImputationBatch.wdl" as Glimpse2SVImputationBatch

workflow Glimpse2SVImputation {
    String pipeline_version = "0.0.1"
    String preprocess_pls_gvcf_pipeline_version = "0.0.1"
    String batch_pipeline_version = "0.0.1"

    input {
        # inputs for Preprocessign wdl
        File? input_gvcfs_fofn
        File? input_gvcf_idxs_fofn
        File? sample_names_file          # order of sample names must match that of gVCFs

        Array[File]? input_gvcfs
        Array[File]? input_gvcf_idxs
        Array[String]? entity_ids
        File? sample_names_map_file           # TSV map of entity_id (research_id) to id2 for AoU DRAGEN gVCFs;
                                              # Terra struggles with id2 as they are parsed as mixed strings/numbers

        String output_prefix

        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        File? extract_bubble_likelihoods_script
        File? extract_bubble_likelihoods_cargo_toml
        File? extract_bubble_likelihoods_binary
        String? extract_bubble_likelihoods_extra_args

        File paste_vcfs_binary
        Array[String] paste_regions

        # inputs for Batch wdl
        File pop_glimpse2_panel_resources_json
        File? pop_glimpse2_script               # heavily modified version of convert-to-biallelic.py
        File? pop_glimpse2_cargo_toml
        File? pop_glimpse2_binary

        String chromosome

        File genetic_maps_tsv
        File chunked_panel_json

        File? remap_sample_names_file    # TSV with old_name new_name mappings
    }

    call PreprocessPLsGVCF.PreprocessPLsGVCF as PreProcessGVCFs {
        input:
        input_gvcfs_fofn = input_gvcfs_fofn,
        input_gvcf_idxs_fofn = input_gvcf_idxs_fofn,
        sample_names_file = sample_names_file,
        input_gvcfs = input_gvcfs,
        input_gvcf_idxs = input_gvcf_idxs,
        entity_ids = entity_ids,
        sample_names_map_file = sample_names_map_file,
        output_prefix = output_prefix,
        preprocess_panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
        preprocess_panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
        extract_bubble_likelihoods_script = extract_bubble_likelihoods_script,
        extract_bubble_likelihoods_cargo_toml = extract_bubble_likelihoods_cargo_toml,
        extract_bubble_likelihoods_binary = extract_bubble_likelihoods_binary,
        extract_bubble_likelihoods_extra_args = extract_bubble_likelihoods_extra_args,
        paste_vcfs_binary = paste_vcfs_binary,
        paste_regions = paste_regions

    }

    call Glimpse2SVImputationBatch.Glimpse2SVImputationBatch {
        input:
            input_preprocessed_joint_vcf = PreProcessGVCFs.preprocessed_pls_vcf,
            input_preprocessed_joint_vcf_idx = PreProcessGVCFs.preprocessed_pls_vcf_idx,
            remap_sample_names_file = remap_sample_names_file,
            chromosome = chromosome,
            genetic_maps_tsv = genetic_maps_tsv,
            chunked_panel_json = chunked_panel_json,
            output_prefix = output_prefix,
            pop_glimpse2_panel_resources_json = pop_glimpse2_panel_resources_json,
            pop_glimpse2_script = pop_glimpse2_script,
            pop_glimpse2_cargo_toml = pop_glimpse2_cargo_toml,
            pop_glimpse2_binary = pop_glimpse2_binary
    }

    output {
        File glimpse2_bubble_posteriors_vcf = Glimpse2SVImputationBatch.glimpse2_bubble_posteriors_vcf
        File glimpse2_bubble_posteriors_vcf_idx = Glimpse2SVImputationBatch.glimpse2_bubble_posteriors_vcf_idx
        File glimpse2_popped_posteriors_vcf = Glimpse2SVImputationBatch.glimpse2_popped_posteriors_vcf
        File glimpse2_popped_posteriors_vcf_idx = Glimpse2SVImputationBatch.glimpse2_popped_posteriors_vcf_idx
    }
}


