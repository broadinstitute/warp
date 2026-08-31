version 1.0

import "./Glimpse2LowPassImputationBatch.wdl" as Glimpse2LowPassImputationBatch
import "../../../../tasks/wdl/Glimpse2LowPassImputationTasks.wdl" as Glimpse2LowPassImputationTasks

workflow Glimpse2LowPassImputation {
    String pipeline_version = "1.1.0"
    String batch_pipeline_version = "1.1.0"
    String quota_consumed_version = "1.0.1"
    String input_qc_version = "1.1.0"

    input {
        # if multiple data types are provided, the workflow will prioritize cram/cram_indices first, then cram manifest
        Array[File]? crams
        Array[File]? cram_indices
        File? cram_manifest
        String output_basename
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

        # Explicit regions for bcftools extraction (Map of contig -> Array of region strings)
        Map[String, Array[String]] bcftools_shard_map

        # Batch sizes for the 2-level hierarchical merge of shards
        Array[Int] hierarchical_merge_batch_sizes = [500, 50]

        # batch size used by this gateway workflow to split very large sample lists
        Int sample_batch_size = 500

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

    scatter(batch_idx in range(length(SplitIntoSampleBatches.cram_manifest_batches))) {
        call Glimpse2LowPassImputationBatch.Glimpse2LowPassImputationBatch as RunBatch {
            input:
                contigs = contigs,
                reference_panel_prefix = reference_panel_prefix,
                cram_manifest = SplitIntoSampleBatches.cram_manifest_batches[batch_idx],
                fasta = fasta,
                fasta_index = fasta_index,
                output_basename = output_basename + ".batch_" + batch_idx,
                impute_reference_only_variants = impute_reference_only_variants,
                call_indels = call_indels,
                bcftools_shard_map = bcftools_shard_map,
                hierarchical_merge_batch_sizes = hierarchical_merge_batch_sizes,
                glimpse_phase_cpu_override = glimpse_phase_cpu_override,
                gatk_docker = gatk_docker,
                glimpse_docker = glimpse_docker
        }
    }
}
