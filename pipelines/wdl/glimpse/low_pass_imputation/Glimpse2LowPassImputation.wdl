version 1.0

import "./Glimpse2LowPassImputationBatch.wdl" as Glimpse2LowPassImputationBatch
import "../../../../tasks/wdl/Glimpse2LowPassImputationTasks.wdl" as Glimpse2LowPassImputationTasks

workflow Glimpse2LowPassImputation {
    String pipeline_version = "0.0.11"
    String batch_pipeline_version = "0.0.3"
    String quota_consumed_version = "0.0.2"
    String input_qc_version = "1.0.1"

    input {
        # if multiple data types are provided, the workflow will prioritize cram_manifest first, then crams/cram_indices/sample_ids
        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? cram_manifest
        String output_basename

        Array[String] contigs
        # this is the path to a directory that contains sites vcf, sites table, and reference chunks file. should end with a "/"
        String reference_panel_prefix
        File fasta
        File fasta_index
        File ref_dict

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false

        # batch size used when calling SplitIntoBatches to make variant calls from the crams
        Int calling_batch_size = 100

        # batch size used by this gateway workflow to split very large sample lists
        Int sample_batch_size = 1000

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        String glimpse_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse@sha256:a0151730cefaaa9ef78b7f9644c63ebb00ce6cd470fa0d60349daa5eee020aec"
        String docker_merge = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
    }

    if (defined(cram_manifest)) {
        call Glimpse2LowPassImputationTasks.ConvertCramManifestToInputArrays {
            input:
                cram_manifest = select_first([cram_manifest])
        }
    }

    # if neither crams (and cram_indices and sample_ids) nor cram_manifest is provided the workflow will fail at runtime
    Array[String] crams_to_use = select_first([ConvertCramManifestToInputArrays.crams, crams])
    Array[String] cram_indices_to_use = select_first([ConvertCramManifestToInputArrays.cram_indices, cram_indices])
    Array[String] sample_ids_to_use = select_first([ConvertCramManifestToInputArrays.sample_ids, sample_ids])

    call Glimpse2LowPassImputationBatch.SplitIntoBatches as SplitIntoSampleBatches {
        input:
            batch_size = sample_batch_size,
            crams = crams_to_use,
            cram_indices = cram_indices_to_use,
            sample_ids = sample_ids_to_use
    }

    scatter(batch_idx in range(length(SplitIntoSampleBatches.crams_batches))) {
        Int batch_sample_count = length(SplitIntoSampleBatches.sample_ids_batches[batch_idx])
        call Glimpse2LowPassImputationBatch.Glimpse2LowPassImputationBatch as RunBatch {
            input:
                contigs = contigs,
                reference_panel_prefix = reference_panel_prefix,
                crams = SplitIntoSampleBatches.crams_batches[batch_idx],
                cram_indices = SplitIntoSampleBatches.cram_indices_batches[batch_idx],
                sample_ids = SplitIntoSampleBatches.sample_ids_batches[batch_idx],
                fasta = fasta,
                fasta_index = fasta_index,
                output_basename = output_basename + ".batch_" + batch_idx,
                ref_dict = ref_dict,
                impute_reference_only_variants = impute_reference_only_variants,
                call_indels = call_indels,
                calling_batch_size = calling_batch_size,
                gatk_docker = gatk_docker,
                glimpse_docker = glimpse_docker
        }
    }

    # Transpose from [batch x contig] to [contig x batch], then paste sample columns together per contig.
    # Use the raw ligated VCFs (all sites) — every batch is imputed to the same reference panel sites,
    # so site lists are identical across batches and the paste md5sum check will pass.
    scatter(contig_idx in range(length(contigs))) {
        Array[File] batch_vcfs_for_contig = transpose(RunBatch.imputed_contig_ligated_vcfs)[contig_idx]
        Array[File] batch_vcf_indices_for_contig = transpose(RunBatch.imputed_contig_ligated_vcf_indices)[contig_idx]

        # For multiple batches: paste sample columns together, recompute AF/INFO as weighted averages.
        # For a single batch: skip both steps since the batch VCF already has correct annotations.
        if (length(SplitIntoSampleBatches.crams_batches) > 1) {
            # Extract AF and INFO annotations from each batch before merging so they can be recalculated
            scatter(batch_annot_idx in range(length(batch_vcfs_for_contig))) {
                call Glimpse2LowPassImputationTasks.ExtractAnnotations {
                    input:
                        imputed_vcf = batch_vcfs_for_contig[batch_annot_idx],
                        imputed_vcf_index = batch_vcf_indices_for_contig[batch_annot_idx],
                        batch_index = batch_annot_idx,
                        docker_extract_annotations = gatk_docker
                }
            }

            call Glimpse2LowPassImputationTasks.MergeSampleChunksVcfsWithPaste as MergeContigVcfs {
                input:
                    input_vcfs = batch_vcfs_for_contig,
                    output_vcf_basename = output_basename + "." + contigs[contig_idx] + ".imputed.merged"
            }

            call Glimpse2LowPassImputationTasks.RecomputeAndAnnotate {
                input:
                    merged_vcf = MergeContigVcfs.output_vcf,
                    annotations = ExtractAnnotations.annotations,
                    num_samples = batch_sample_count,
                    output_basename = output_basename + "." + contigs[contig_idx] + ".imputed.merged.reannotated",
                    docker_merge = docker_merge
            }
        }

        File annotated_contig_vcf = select_first([RecomputeAndAnnotate.merged_imputed_vcf, batch_vcfs_for_contig[0]])

        # Now that the full cohort is merged and annotations are correct, split into variant-only and hom-ref-only
        call Glimpse2LowPassImputationTasks.SelectVariantRecordsOnly as SelectContigVariants {
            input:
                vcf = annotated_contig_vcf,
                basename = output_basename + "." + contigs[contig_idx] + ".imputed.merged.only_variants"
        }

        call Glimpse2LowPassImputationTasks.CreateHomRefSitesOnlyVcf as CreateContigHomRefVcf {
            input:
                vcf = annotated_contig_vcf,
                basename = output_basename + "." + contigs[contig_idx] + ".imputed.merged.only_hom_ref.sites_only"
        }
    }

    Array[File] contig_variant_vcfs = SelectContigVariants.output_vcf
    Array[File] contig_hom_ref_vcfs = CreateContigHomRefVcf.output_vcf
    Array[File] batch_coverage_metrics = select_all(RunBatch.coverage_metrics)

    if (length(batch_coverage_metrics) > 0) {
        call Glimpse2LowPassImputationTasks.MergeCoverageMetrics as MergeBatchCoverageMetrics {
            input:
                coverage_metrics = batch_coverage_metrics,
                output_basename = output_basename
        }
    }

    call Glimpse2LowPassImputationBatch.GatherVcfsNoIndex {
        input:
            input_vcfs = contig_variant_vcfs,
            output_vcf_basename = output_basename + ".imputed",
            gatk_docker = gatk_docker
    }

    call Glimpse2LowPassImputationBatch.CreateVcfIndexAndMd5 {
        input:
            vcf_input = GatherVcfsNoIndex.output_vcf,
            gatk_docker = gatk_docker,
            preemptible = 0
    }

    call Glimpse2LowPassImputationBatch.GatherVcfsNoIndex as GatherVcfsNoIndexHomRefOnly {
        input:
            input_vcfs = contig_hom_ref_vcfs,
            output_vcf_basename = output_basename + ".imputed.hom_ref_sites_only",
            gatk_docker = gatk_docker
    }

    call Glimpse2LowPassImputationBatch.CreateVcfIndexAndMd5 as CreateVcfIndexAndMd5HomRefOnly {
        input:
            vcf_input = GatherVcfsNoIndexHomRefOnly.output_vcf,
            gatk_docker = gatk_docker,
            preemptible = 0
    }

    call Glimpse2LowPassImputationBatch.CollectQCMetrics {
        input:
            imputed_vcf = GatherVcfsNoIndex.output_vcf,
            output_basename = output_basename
    }

    output {
        File imputed_vcf = CreateVcfIndexAndMd5.output_vcf
        File imputed_vcf_index = CreateVcfIndexAndMd5.output_vcf_index
        File imputed_vcf_md5sum = CreateVcfIndexAndMd5.output_vcf_md5sum

        File imputed_hom_ref_sites_only_vcf = CreateVcfIndexAndMd5HomRefOnly.output_vcf
        File imputed_hom_ref_sites_only_vcf_index = CreateVcfIndexAndMd5HomRefOnly.output_vcf_index
        File imputed_hom_ref_sites_only_vcf_md5 = CreateVcfIndexAndMd5HomRefOnly.output_vcf_md5sum

        File qc_metrics = CollectQCMetrics.qc_metrics
        File? coverage_metrics = MergeBatchCoverageMetrics.merged_coverage_metrics
    }
}
