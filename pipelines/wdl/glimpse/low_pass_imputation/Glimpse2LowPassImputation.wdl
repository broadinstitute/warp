version 1.0

import "./Glimpse2LowPassImputationBatch.wdl" as Glimpse2LowPassImputationBatch
import "../../../../tasks/wdl/Glimpse2LowPassImputationTasks.wdl" as Glimpse2LowPassImputationTasks

workflow Glimpse2LowPassImputation {
    String pipeline_version = "2.1.2"
    String batch_pipeline_version = "1.1.2"
    String quota_consumed_version = "1.0.2"
    String input_qc_version = "1.1.2"

    input {
        # if multiple data types are provided, the workflow will prioritize cram/cram_indices first, then cram manifest
        Array[File]? crams
        Array[File]? cram_indices
        File? cram_manifest
        String output_basename

        # remove RAF when RAF <= raf_censoring_threshold
        Float raf_censoring_threshold = 0.0

        # Optional filter: variants with INFO score below this threshold will be excluded from the final output VCF
        Float info_filter_for_inclusion = 0.0

        Array[String] contigs
        # this is the path to a directory that contains sites vcf, sites table, and reference chunks file. should end with a "/"
        String reference_panel_prefix
        File fasta
        File fasta_index
        File ref_dict

        # optional additional header line to add to the output VCF
        String? pipeline_header_line

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false

        # batch size used when calling SplitIntoBatches to make variant calls from the crams
        Int calling_batch_size = 100

        # batch size used by this gateway workflow to split very large sample lists
        Int sample_batch_size = 500

        # size of postprocessing shards
        Int postprocessing_shard_size_bp = 25000000

        # override for cpu used for glimpse phase task. Mostly used to set to 1 for determinism in testing, defaults to 4
        Int? glimpse_phase_cpu_override

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        String glimpse_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.3.0-8671138-1785933808"
        String docker_merge = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
    }

    if (defined(crams)) {
        call Glimpse2LowPassImputationTasks.ConvertInputArraysToManifest {
            input:
                cram_paths = select_first([crams]),
                cram_index_paths = select_first([cram_indices])
        }
    }

    # if neither crams (and cram_indices) nor cram_manifest is provided the workflow will fail at runtime
    File cram_manifest_to_use = select_first([ConvertInputArraysToManifest.output_manifest, cram_manifest])

    call Glimpse2LowPassImputationTasks.SplitCramManifestIntoBatches as SplitIntoSampleBatches {
        input:
            batch_size = sample_batch_size,
            cram_manifest = cram_manifest_to_use
    }
    
    Boolean multiple_batches = length(SplitIntoSampleBatches.cram_manifest_batches) > 1

    scatter(batch_idx in range(length(SplitIntoSampleBatches.cram_manifest_batches))) {
        call Glimpse2LowPassImputationBatch.Glimpse2LowPassImputationBatch as RunBatch {
            input:
                contigs = contigs,
                reference_panel_prefix = reference_panel_prefix,
                cram_manifest = SplitIntoSampleBatches.cram_manifest_batches[batch_idx],
                fasta = fasta,
                fasta_index = fasta_index,
                output_basename = output_basename + ".batch_" + batch_idx,
                ref_dict = ref_dict,
                pipeline_header_line = pipeline_header_line,
                impute_reference_only_variants = impute_reference_only_variants,
                call_indels = call_indels,
                calling_batch_size = calling_batch_size,
                glimpse_phase_cpu_override = glimpse_phase_cpu_override,
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

        # For multiple batches: extract annotations, paste sample columns together, then recompute AF/INFO as weighted averages and apply filters.
        if (multiple_batches) {
            call Glimpse2LowPassImputationTasks.CreateContigRegions {
                input:
                    fasta_index = fasta_index,
                    contig = contigs[contig_idx],
                    shard_size_bp = postprocessing_shard_size_bp
            }

            scatter(region in CreateContigRegions.regions) {
                # Format region for safe filenames
                String region_id = sub(region, ":", "_")
                
                # Extract AF and INFO annotations from each batch before merging so they can be recalculated
                scatter(batch_annot_idx in range(length(batch_vcfs_for_contig))) {
                    call Glimpse2LowPassImputationTasks.ExtractAnnotations {
                        input:
                            imputed_vcf = batch_vcfs_for_contig[batch_annot_idx],
                            imputed_vcf_index = batch_vcf_indices_for_contig[batch_annot_idx],
                            batch_index = batch_annot_idx,
                            region = region,
                            docker_extract_annotations = docker_merge
                    }
                }

                call Glimpse2LowPassImputationTasks.MergeSampleChunksVcfsWithPaste as MergeContigVcfs {
                    input:
                        input_vcfs = batch_vcfs_for_contig,
                        input_vcf_indices = batch_vcf_indices_for_contig,
                        region = region,
                        output_vcf_basename = output_basename + "." + contigs[contig_idx] + "." + region_id + ".imputed.merged"
                }

                call Glimpse2LowPassImputationTasks.RecomputeAndAnnotate {
                    input:
                        merged_vcf = MergeContigVcfs.output_vcf,
                        merged_vcf_index = MergeContigVcfs.output_vcf_index,
                        annotations = ExtractAnnotations.annotations,
                        num_samples = RunBatch.total_samples,
                        info_filter_threshold = info_filter_for_inclusion,
                        raf_censoring_threshold = raf_censoring_threshold,
                        region = region,
                        output_basename = output_basename + "." + contigs[contig_idx] + "." + region_id + ".imputed.merged.reannotated",
                        docker_merge = docker_merge
                }
            }

            call Glimpse2LowPassImputationTasks.ConcatAndFinalizeVcfs as ConcatAndFinalizeMultiBatch {
                input:
                    vcfs = RecomputeAndAnnotate.merged_imputed_vcf,
                    output_basename = output_basename + "." + contigs[contig_idx] + ".imputed",
                    extra_args = "--naive"
            }
        }

        # For a single batch: apply info filters and RAF censoring via bcftools pipe.
        if (!multiple_batches) {
            call Glimpse2LowPassImputationTasks.FilterInfoAndCensorRafAndFinalizeVcf {
                input:
                    vcf = batch_vcfs_for_contig[0],
                    info_filter_threshold = info_filter_for_inclusion,
                    raf_censoring_threshold = raf_censoring_threshold,
                    output_basename = output_basename + "." + contigs[contig_idx] + ".imputed"
            }
        }

        File finalized_contig_vcf = select_first([ConcatAndFinalizeMultiBatch.concatenated_vcf, FilterInfoAndCensorRafAndFinalizeVcf.filtered_vcf])
        File finalized_contig_vcf_idx = select_first([ConcatAndFinalizeMultiBatch.concatenated_vcf_idx, FilterInfoAndCensorRafAndFinalizeVcf.filtered_vcf_idx])
        File finalized_contig_vcf_md5sum = select_first([ConcatAndFinalizeMultiBatch.concatenated_vcf_md5sum, FilterInfoAndCensorRafAndFinalizeVcf.filtered_vcf_md5sum])
    }

    call Glimpse2LowPassImputationTasks.GatherVcfsNoIndex as GatherContigVcfsForQcMetrics {
        input:
            input_vcfs = finalized_contig_vcf,
            output_vcf_basename = output_basename + ".allchr.imputed",
            gatk_docker = gatk_docker
    }

    call Glimpse2LowPassImputationTasks.CollectQCMetrics {
        input:
            imputed_vcf = GatherContigVcfsForQcMetrics.output_vcf,
            output_basename = output_basename
    }

    output {
        Array[File] imputed_vcfs = finalized_contig_vcf
        Array[File] imputed_vcf_indexes = finalized_contig_vcf_idx
        Array[File] imputed_vcf_md5sums = finalized_contig_vcf_md5sum

        File qc_metrics = CollectQCMetrics.qc_metrics
    }
}
