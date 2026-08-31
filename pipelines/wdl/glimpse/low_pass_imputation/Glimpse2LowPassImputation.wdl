version 1.0

import "./Glimpse2LowPassImputationBatch.wdl" as Glimpse2LowPassImputationBatch
import "../../../../tasks/wdl/Glimpse2LowPassImputationTasks.wdl" as Glimpse2LowPassImputationTasks

workflow Glimpse2LowPassImputation {
    String pipeline_version = "1.1.0"
    String batch_pipeline_version = "1.1.0"
    String quota_consumed_version = "1.0.1"
    String input_qc_version = "1.1.0"

    input {
        Array[File]? crams
        Array[File]? cram_indices
        File? cram_manifest
        String output_basename
        Float info_filter_for_inclusion = 0.0

        Array[String] contigs
        String reference_panel_prefix
        File fasta
        File fasta_index
        File ref_dict

        String? pipeline_header_line

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false

        Map[String, Array[String]] bcftools_shard_map
        Int hierarchical_merge_batch_size = 50

        Int sample_batch_size = 500
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
                hierarchical_merge_batch_size = hierarchical_merge_batch_size,
                glimpse_phase_cpu_override = glimpse_phase_cpu_override,
                gatk_docker = gatk_docker,
                glimpse_docker = glimpse_docker
        }
    }
}
