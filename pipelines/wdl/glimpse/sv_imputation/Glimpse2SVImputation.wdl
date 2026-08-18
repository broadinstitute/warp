version 1.0

import "./PreprocessPLsGVCF.wdl" as PreprocessPLsGVCF
import "./Glimpse2SVImputationBatch.wdl" as Glimpse2SVImputationBatch
import "../../../../tasks/wdl/Glimpse2SVImputationTasks.wdl" as Glimpse2SVImputationTasks

workflow Glimpse2SVImputation {
    String pipeline_version = "0.0.21"
    String preprocess_pls_gvcf_pipeline_version = "0.0.12"
    String batch_pipeline_version = "0.0.16"

    input {
        # if both array inputs and gvcf_manifest are provided, array inputs take precedence
        Array[File]? input_gvcfs
        Array[File]? input_gvcf_idxs
        Array[String]? sample_ids
        File? gvcf_manifest
        Int sample_batch_size = 1000

        String output_prefix

        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        String? extract_bubble_likelihoods_extra_args

        Array[String] paste_regions

        # inputs for Batch wdl
        Array[String] chromosomes
        File genetic_maps_tsv
        File ref_dict
        File chunked_panel_json

        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites --Kpbwt 1000 --main 10 --burnin 5 --err-imp 1E-3"

        # override for cpu used for glimpse phase task. Mostly used to set to 1 for determinism in testing, defaults to 4
        Int? glimpse_phase_cpu_override

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json

        # Optional filter: variants with INFO score below this threshold will be excluded from the final output VCFs
        Float info_filter_for_inclusion = 0.0

        # optional additional header line to add to the output VCF
        String? pipeline_header_line

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.2.0-8671138-1784681771"
        String merge_docker = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }

    Boolean using_arrays = defined(input_gvcfs) && defined(input_gvcf_idxs) && defined(sample_ids)

    if (using_arrays) {
        call Glimpse2SVImputationTasks.ConvertInputArraysToManifest {
            input:
                gvcf_paths = select_first([input_gvcfs]),
                gvcf_index_paths = select_first([input_gvcf_idxs]),
                sample_ids = select_first([sample_ids])
        }
    }

    # if neither the full array input set nor gvcf_manifest is provided the workflow will fail at runtime
    File gvcf_manifest_to_use = select_first([ConvertInputArraysToManifest.output_manifest, gvcf_manifest])

    call Glimpse2SVImputationTasks.SplitVcfManifestIntoBatches as SplitIntoSampleBatches {
        input:
            batch_size = sample_batch_size,
            gvcf_manifest = gvcf_manifest_to_use
    }

    scatter (batch_idx in range(length(SplitIntoSampleBatches.gvcf_manifest_batches))) {
        call PreprocessPLsGVCF.PreprocessPLsGVCF as PreProcessGVCFsBatch {
            input:
                input_gvcf_manifest = SplitIntoSampleBatches.gvcf_manifest_batches[batch_idx],
                preprocess_panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
                preprocess_panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
                extract_bubble_likelihoods_extra_args = extract_bubble_likelihoods_extra_args,
                paste_regions = paste_regions
        }

        call Glimpse2SVImputationBatch.Glimpse2SVImputationBatch as RunBatch {
            input:
                input_preprocessed_joint_vcf = PreProcessGVCFsBatch.preprocessed_pls_vcf,
                input_preprocessed_joint_vcf_idx = PreProcessGVCFsBatch.preprocessed_pls_vcf_idx,
                chromosomes = chromosomes,
                genetic_maps_tsv = genetic_maps_tsv,
                ref_dict = ref_dict,
                chunked_panel_json = chunked_panel_json,
                extra_phase_args = extra_phase_args,
                output_prefix = output_prefix + ".batch_" + batch_idx,
                pop_glimpse2_panel_resources_json = pop_glimpse2_panel_resources_json,
                glimpse2_docker = glimpse2_docker,
                glimpse_phase_cpu_override = glimpse_phase_cpu_override,
                pipeline_header_line = pipeline_header_line
        }
    }

    scatter (contig_idx in range(length(chromosomes))) {
        Array[File] popped_bcfs_for_contig = transpose(RunBatch.glimpse2_popped_posteriors_vcf)[contig_idx]
        Array[File] popped_bcf_idxs_for_contig = transpose(RunBatch.glimpse2_popped_posteriors_vcf_idx)[contig_idx]

        if (length(SplitIntoSampleBatches.gvcf_manifest_batches) > 1) {
            scatter (batch_annot_idx in range(length(popped_bcfs_for_contig))) {
                call Glimpse2SVImputationTasks.ExtractAnnotations as ExtractPoppedAnnotations {
                    input:
                        imputed_vcf = popped_bcfs_for_contig[batch_annot_idx],
                        imputed_vcf_index = popped_bcf_idxs_for_contig[batch_annot_idx],
                        batch_index = batch_annot_idx,
                        docker_extract_annotations = gatk_docker
                }
            }

            call Glimpse2SVImputationTasks.MergeSampleChunksVcfsWithPaste as MergePoppedContigVcfs {
                input:
                    input_vcfs = popped_bcfs_for_contig,
                    output_vcf_basename = output_prefix + "." + chromosomes[contig_idx] + ".glimpse2.popped.merged"
            }

            call Glimpse2SVImputationTasks.RecomputeAndAnnotate as RecomputePoppedAfInfo {
                input:
                    merged_vcf = MergePoppedContigVcfs.output_vcf,
                    annotations = ExtractPoppedAnnotations.annotations,
                    num_samples = PreProcessGVCFsBatch.num_samples,
                    output_basename = output_prefix + "." + chromosomes[contig_idx] + ".glimpse2.popped.merged.reannotated",
                    docker_merge = merge_docker
            }
        }

        File final_popped_contig_vcf = select_first([RecomputePoppedAfInfo.merged_imputed_vcf, popped_bcfs_for_contig[0]])

        if (info_filter_for_inclusion > 0.0) {
            call Glimpse2SVImputationTasks.FilterVcfByInfo as FilterPoppedContigByInfo {
                input:
                    vcf = final_popped_contig_vcf,
                    info_threshold = info_filter_for_inclusion,
                    output_prefix = output_prefix + "." + chromosomes[contig_idx] + ".glimpse2.popped.info_filtered"
            }
        }

        File final_filtered_popped_contig_vcf = select_first([FilterPoppedContigByInfo.output_vcf, final_popped_contig_vcf])

        call Glimpse2SVImputationTasks.CreateVcfIndexAndMd5 as IndexFinalPoppedContig {
            input:
                vcf_input = final_filtered_popped_contig_vcf,
                output_basename = output_prefix + "." + chromosomes[contig_idx],
                gatk_docker = gatk_docker,
                preemptible = 0
        }
    }

    output {
        Array[File] glimpse2_popped_posteriors_vcf = IndexFinalPoppedContig.output_vcf
        Array[File] glimpse2_popped_posteriors_vcf_idx = IndexFinalPoppedContig.output_vcf_index
    }
}

