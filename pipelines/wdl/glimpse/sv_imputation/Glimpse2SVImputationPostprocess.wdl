version 1.0

import "./PreprocessPLsGVCF.wdl" as PreprocessPLsGVCF
import "./Glimpse2SVImputationPostprocessBatch.wdl" as Glimpse2SVImputationBatch
import "./MultilevelHierarchicallyPasteVcfsStreaming.wdl" as MultilevelMerge
import "../../../../tasks/wdl/Glimpse2SVImputationTasks.wdl" as Glimpse2SVImputationTasks

workflow Glimpse2SVImputation {
    String pipeline_version = "0.0.28"
    String preprocess_pls_gvcf_pipeline_version = "0.0.16"
    String batch_pipeline_version = "0.0.20"
    String quota_consumed_version = "0.0.2"
    String input_qc_version = "0.0.2"

    input {
        # if both array inputs and gvcf_manifest are provided, array inputs take precedence
        Array[File]? input_gvcfs
        Array[File]? input_gvcf_idxs
        File? gvcf_manifest
        Int sample_batch_size = 1000

        String output_basename

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

        # Optional overrides to skip batch-level imputation step and jump straight to marginalize collisions
        Array[Array[File]]? batch_chromosome_posteriors_vcfs_or_bcfs
        Array[Array[File]]? batch_chromosome_posteriors_vcf_or_bcf_idxs
        Array[Int]? override_batch_num_samples

        # Multilevel Paste Configuration
        Array[Int] merge_batch_sizes = [100]        # adjust or add levels if needed
        Array[Boolean] merge_do_localization = [true]
        Array[Int] merge_timeouts_min = [120]

        # Optional filter: variants with INFO score below this threshold will be excluded from the final output VCFs
        Float info_filter_for_inclusion = 0.0

        # optional additional header line to add to the output VCF
        String? pipeline_header_line

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.2.0-8671138-1784681771"
        String merge_docker = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }

    Boolean using_arrays = defined(input_gvcfs) && defined(input_gvcf_idxs)

    if (using_arrays) {
        call Glimpse2SVImputationTasks.ConvertInputArraysToManifest {
            input:
                gvcf_paths = select_first([input_gvcfs]),
                gvcf_index_paths = select_first([input_gvcf_idxs])
        }
    }

    # if neither the full array input set nor gvcf_manifest is provided the workflow will fail at runtime
    File gvcf_manifest_to_use = select_first([ConvertInputArraysToManifest.output_gvcf_manifest, gvcf_manifest])

    call Glimpse2SVImputationTasks.SplitVcfManifestIntoBatches as SplitIntoSampleBatches {
        input:
            batch_size = sample_batch_size,
            gvcf_manifest = gvcf_manifest_to_use
    }

    Map[String, ChunkedPanelChromosome] chunked_panel = read_json(chunked_panel_json)
    Map[String, PopAndMarginalizePanelResourcesChromosome] pop_glimpse2_panel_resources = read_json(pop_glimpse2_panel_resources_json)

    scatter (batch_idx in range(length(SplitIntoSampleBatches.gvcf_manifest_batches))) {
        
        Boolean has_posteriors = defined(batch_chromosome_posteriors_vcfs_or_bcfs)

        if (!has_posteriors) {
            call PreprocessPLsGVCF.PreprocessPLsGVCF as PreProcessGVCFsBatch {
                input:
                    input_gvcf_manifest = SplitIntoSampleBatches.gvcf_manifest_batches[batch_idx],
                    preprocess_panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
                    preprocess_panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
                    extract_bubble_likelihoods_extra_args = extract_bubble_likelihoods_extra_args,
                    paste_regions = paste_regions
            }
        }

        Int current_batch_num_samples = if defined(override_batch_num_samples) then select_first([override_batch_num_samples])[batch_idx] else select_first([PreProcessGVCFsBatch.num_samples])

        if (has_posteriors) {
            Array[File] extracted_vcfs = select_first([batch_chromosome_posteriors_vcfs_or_bcfs])[batch_idx]
            Array[File] extracted_idxs = select_first([batch_chromosome_posteriors_vcf_or_bcf_idxs])[batch_idx]
        }

        call Glimpse2SVImputationBatch.Glimpse2SVImputationBatch as RunBatch {
            input:
                input_preprocessed_joint_vcf_or_bcf = select_first([PreProcessGVCFsBatch.preprocessed_pls_bcf, "dummy.bcf"]),
                input_preprocessed_joint_vcf_or_bcf_idx = select_first([PreProcessGVCFsBatch.preprocessed_pls_bcf_idx, "dummy.bcf.csi"]),
                chromosomes = chromosomes,
                genetic_maps_tsv = genetic_maps_tsv,
                ref_dict = ref_dict,
                chunked_panel_json = chunked_panel_json,
                extra_phase_args = extra_phase_args,
                output_basename = output_basename + ".batch_" + batch_idx,
                pop_glimpse2_panel_resources_json = pop_glimpse2_panel_resources_json,
                glimpse2_docker = glimpse2_docker,
                glimpse_phase_cpu_override = glimpse_phase_cpu_override,
                pipeline_header_line = pipeline_header_line,
                chromosome_posteriors_vcfs_or_bcfs = extracted_vcfs,
                chromosome_posteriors_vcf_or_bcf_idxs = extracted_idxs
        }
    }

    scatter (contig_idx in range(length(chromosomes))) {
        String chr = chromosomes[contig_idx]
        Array[File] popped_bcfs_for_contig = transpose(RunBatch.glimpse2_popped_posteriors_vcf)[contig_idx]
        Array[File] popped_bcf_idxs_for_contig = transpose(RunBatch.glimpse2_popped_posteriors_vcf_idx)[contig_idx]

        Array[String] contig_regions = select_first([pop_glimpse2_panel_resources[chr].pop_regions, chunked_panel[chr].output_regions])

        if (length(SplitIntoSampleBatches.gvcf_manifest_batches) > 1) {
            
            scatter (region in contig_regions) {
                scatter (batch_annot_idx in range(length(popped_bcfs_for_contig))) {
                    call Glimpse2SVImputationTasks.ExtractAnnotations as ExtractPoppedAnnotations {
                        input:
                            imputed_vcf_or_bcf = popped_bcfs_for_contig[batch_annot_idx],
                            imputed_vcf_or_bcf_index = popped_bcf_idxs_for_contig[batch_annot_idx],
                            batch_index = batch_annot_idx,
                            region = region,
                            docker_extract_annotations = gatk_docker
                    }
                }

                call MultilevelMerge.MultilevelHierarchicallyMergeVcfs as MergePoppedRegion {
                    input:
                        vcfs_or_bcfs_array = popped_bcfs_for_contig,
                        vcf_or_bcf_idxs_array = popped_bcf_idxs_for_contig,
                        regions = [region],
                        batch_sizes = merge_batch_sizes,
                        do_localization = merge_do_localization,
                        timeouts_min = merge_timeouts_min,
                        output_basename = output_basename + "." + chr + "." + region + ".glimpse2.popped.merged"
                }

                call Glimpse2SVImputationTasks.RecomputeAndAnnotate as RecomputePoppedAfInfo {
                    input:
                        merged_vcf_or_bcf = MergePoppedRegion.merged_bcf,
                        annotations = ExtractPoppedAnnotations.annotations,
                        num_samples = current_batch_num_samples,
                        output_basename = output_basename + "." + chr + "." + region + ".glimpse2.popped.merged.reannotated",
                        docker_merge = merge_docker
                }
            }

            call Glimpse2SVImputationTasks.ConcatBcfs as ConcatContigRegions {
                input:
                    bcfs = RecomputePoppedAfInfo.merged_imputed_bcf,
                    bcf_idxs = RecomputePoppedAfInfo.merged_imputed_bcf_idx,
                    output_basename = output_basename + "." + chr + ".glimpse2.popped.merged.reannotated",
                    extra_args = "--naive"
            }
        }

        File final_popped_contig_vcf = select_first([ConcatContigRegions.concatenated_bcf, popped_bcfs_for_contig[0]])

        call Glimpse2SVImputationTasks.CreateVcfIndexAndMd5 as IndexFinalPoppedContig {
            input:
                vcf_input_or_bcf = final_popped_contig_vcf,
                output_basename = output_basename + "." + chr,
                info_filter_threshold = info_filter_for_inclusion,
                gatk_docker = gatk_docker,
                preemptible = 0
        }
    }

    output {
        Array[File] imputed_vcfs = IndexFinalPoppedContig.output_vcf
        Array[File] imputed_vcf_indexes = IndexFinalPoppedContig.output_vcf_index
    }
}

struct ChunkedPanelChromosome {
    Array[String] input_regions
    Array[String] output_regions
    Array[String] panel_split_chunk_bins
    Array[Int]? phase_base_mem
}

struct PopAndMarginalizePanelResourcesChromosome {
    String panel_bubble_split_sites_only_vcf
    String panel_bubble_split_sites_only_vcf_idx
    String panel_id_split_vcf_gz
    String panel_id_split_vcf_gz_tbi
    Array[String]? pop_regions
}
