version 1.0

# BroadInternalImputation is now deprecated 2025-03-06

import "../../../../../pipelines/wdl/arrays/imputation/Imputation.wdl" as ImputationPipeline
import "../../../../../tasks/wdl/InternalImputationTasks.wdl" as InternalImputationTasks
import "../../../../../tasks/wdl/InternalTasks.wdl" as InternalTasks

workflow BroadInternalImputation {
    meta {
        description: "Push outputs of Imputation.wdl to TDR dataset table ImputationOutputsTable and split out Imputation arrays into ImputationWideOutputsTable."
        allowNestedInputs: true
    }
    String pipeline_version = "1.1.16"
    
    input {
        # inputs to wrapper task 
        String tdr_dataset_id
        String tdr_target_table_name
        String prs_cf_trigger_bucket_path

        # required inputs to Imputation.wdl
        Array[String]   contigs
        File            genetic_maps_eagle
        String          output_callset_name
        File            ref_dict
        String          reference_panel_path

        # optional inputs to Imputation.wdl but required for eMerge
        Array[File]     single_sample_vcfs
        Array[File]     single_sample_vcf_indices
        Array[String]   chip_well_barcodes

        Array[String]   lab_batches
        String          timestamp
    }

    call ImputationPipeline.Imputation {
        input:
            contigs = contigs,
            genetic_maps_eagle          = genetic_maps_eagle,
            output_callset_name         = output_callset_name,
            ref_dict                    = ref_dict,
            reference_panel_path        = reference_panel_path,
            single_sample_vcfs          = single_sample_vcfs,
            single_sample_vcf_indices   = single_sample_vcf_indices
    }

    call InternalImputationTasks.FormatImputationOutputs {
        input:
            imputed_single_sample_vcfs          = Imputation.imputed_single_sample_vcfs,
            imputed_single_sample_vcf_indices   = Imputation.imputed_single_sample_vcf_indices,
            chip_well_barcodes                  = chip_well_barcodes,
            imputed_multisample_vcf             = Imputation.imputed_multisample_vcf,
            imputed_multisample_vcf_index       = Imputation.imputed_multisample_vcf_index,
            aggregated_imputation_metrics       = Imputation.aggregated_imputation_metrics,
            chunks_info                         = Imputation.chunks_info,
            failed_chunks                       = Imputation.failed_chunks,
            n_failed_chunks                     = Imputation.n_failed_chunks
    }

    call InternalTasks.IngestOutputsToTDR as IngestToImputationOutputsTable {
        input:
            tdr_dataset_id          = tdr_dataset_id,
            tdr_target_table_name   = tdr_target_table_name,
            outputs_tsv             = FormatImputationOutputs.ingest_outputs_tsv
    }

    call InternalImputationTasks.FormatImputationWideOutputs {
        input:
            imputed_single_sample_vcfs          = Imputation.imputed_single_sample_vcfs,
            imputed_single_sample_vcf_indices   = Imputation.imputed_single_sample_vcf_indices
    }

    call InternalTasks.IngestOutputsToTDR as IngestToImputationWideOutputsTable {
        input:
            tdr_dataset_id          = tdr_dataset_id,
            tdr_target_table_name   = "ImputationWideOutputsTable",
            outputs_tsv             = FormatImputationWideOutputs.ingest_outputs_wide_tsv,
            in_load_tag             = IngestToImputationOutputsTable.load_tag
    }

    call InternalImputationTasks.TriggerPrsWithImputationTsv {
        input:
            run_task                = IngestToImputationWideOutputsTable.ingest_logs,
            imputation_outputs_tsv  = FormatImputationOutputs.ingest_outputs_tsv,
            trigger_bucket_path     = prs_cf_trigger_bucket_path,
            timestamp               = timestamp,
            lab_batches             = lab_batches
    }

    output {
        File aggregated_imputation_metrics              = Imputation.aggregated_imputation_metrics
        File chunks_info                                = Imputation.chunks_info
        File failed_chunks                              = Imputation.failed_chunks
        File imputed_multisample_vcf                    = Imputation.imputed_multisample_vcf
        File imputed_multisample_vcf_index              = Imputation.imputed_multisample_vcf_index
        Array[File]? imputed_single_sample_vcfs         = Imputation.imputed_single_sample_vcfs
        Array[File]? imputed_single_sample_vcf_indices  = Imputation.imputed_single_sample_vcf_indices
        File n_failed_chunks                            = Imputation.n_failed_chunks
    }
}
