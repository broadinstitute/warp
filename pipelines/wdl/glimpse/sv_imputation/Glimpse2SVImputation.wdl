version 1.0

import "./PreprocessPLsGVCF.wdl" as PreprocessPLsGVCF
import "./Glimpse2SVImputationBatch.wdl" as Glimpse2SVImputationBatch

workflow Glimpse2SVImputation {
    String pipeline_version = "0.0.10"
    String preprocess_pls_gvcf_pipeline_version = "0.0.6"
    String batch_pipeline_version = "0.0.8"

    input {
        # inputs for Preprocessign wdl
        File? input_gvcfs_fofn
        File? input_gvcf_idxs_fofn
        File? sample_ids_file          # order of sample ids must match that of gVCFs

        Array[File]? input_gvcfs
        Array[File]? input_gvcf_idxs
        Array[String]? sample_ids

        String output_prefix

        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        String? extract_bubble_likelihoods_extra_args

        Array[String] paste_regions

        # inputs for Batch wdl
        Array[String] chromosomes
        File genetic_maps_tsv
        File chunked_panel_json

        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites --Kpbwt 1000 --main 10 --burnin 5 --err-imp 1E-3"
        
        # override for cpu used for glimpse phase task. Mostly used to set to 1 for determinism in testing, defaults to 4
        Int? glimpse_phase_cpu_override

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.2.0-8671138-1784681771"
    }

    call PreprocessPLsGVCF.PreprocessPLsGVCF as PreProcessGVCFs {
        input:
        input_gvcfs_fofn = input_gvcfs_fofn,
        input_gvcf_idxs_fofn = input_gvcf_idxs_fofn,
        sample_ids_file = sample_ids_file,
        input_gvcfs = input_gvcfs,
        input_gvcf_idxs = input_gvcf_idxs,
        sample_ids = sample_ids,
        preprocess_panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
        preprocess_panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
        extract_bubble_likelihoods_extra_args = extract_bubble_likelihoods_extra_args,
        paste_regions = paste_regions

    }

    call Glimpse2SVImputationBatch.Glimpse2SVImputationBatch {
        input:
            input_preprocessed_joint_vcf = PreProcessGVCFs.preprocessed_pls_vcf,
            input_preprocessed_joint_vcf_idx = PreProcessGVCFs.preprocessed_pls_vcf_idx,
            chromosomes = chromosomes,
            genetic_maps_tsv = genetic_maps_tsv,
            chunked_panel_json = chunked_panel_json,
            extra_phase_args = extra_phase_args,
            output_prefix = output_prefix,
            pop_glimpse2_panel_resources_json = pop_glimpse2_panel_resources_json,
            glimpse2_docker = glimpse2_docker,
            glimpse_phase_cpu_override = glimpse_phase_cpu_override
    }

    output {
        Array[File] glimpse2_bubble_posteriors_vcf = Glimpse2SVImputationBatch.glimpse2_bubble_posteriors_vcf
        Array[File] glimpse2_bubble_posteriors_vcf_idx = Glimpse2SVImputationBatch.glimpse2_bubble_posteriors_vcf_idx
        Array[File] glimpse2_popped_posteriors_vcf = Glimpse2SVImputationBatch.glimpse2_popped_posteriors_vcf
        Array[File] glimpse2_popped_posteriors_vcf_idx = Glimpse2SVImputationBatch.glimpse2_popped_posteriors_vcf_idx
    }
}


