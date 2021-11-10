version 1.0

import "Imputation.wdl" as ImputationPipeline

workflow imputation_outputs_to_TDR {
    meta {
        description: "Push outputs of Imputation.wdl to TDR dataset table ImputationOutputsTable."
    }

    input {
        # inputs to wrapper task 
        String workspace_name
        String workspace_bucket
        String tdr_dataset_id
        String tdr_target_table_name

        # required inputs to Imputation.wdl
        Array[String]   contigs
        File            genetic_maps_eagle
        File            haplotype_database
        String          output_callset_name
        File            ref_dict
        String          reference_panel_path

        # optional inputs to Imputation.wdl but required for eMerge
        Array[File]     single_sample_vcfs
        Array[File]     single_sample_vcf_indices
    }

    call ImputationPipeline.Imputation {
        input:
            contigs = contigs,
            genetic_maps_eagle          = genetic_maps_eagle,
            haplotype_database          = haplotype_database,
            output_callset_name         = output_callset_name,
            ref_dict                    = ref_dict,
            reference_panel_path        = reference_panel_path,
            single_sample_vcfs          = single_sample_vcfs,
            single_sample_vcf_indices   = single_sample_vcf_indices
    }

    call format_imputation_outputs {
        input:
            imputed_single_sample_vcfs          = Imputation.imputed_single_sample_vcfs,
            imputed_single_sample_vcf_indices   = Imputation.imputed_single_sample_vcf_indices,
            imputed_multisample_vcf             = Imputation.imputed_multisample_vcf,
            imputed_multisample_vcf_index       = Imputation.imputed_multisample_vcf_index,
            aggregated_imputation_metrics       = Imputation.aggregated_imputation_metrics,
            chunks_info                         = Imputation.chunks_info,
            failed_chunks                       = Imputation.failed_chunks,
            n_failed_chunks                     = Imputation.n_failed_chunks
    }

    call ingest_outputs_to_tdr {
        input:
            workspace_name          = workspace_name,
            workspace_bucket        = workspace_bucket,
            tdr_dataset_id          = tdr_dataset_id,
            tdr_target_table_name   = tdr_target_table_name,
            outputs_tsv             = format_imputation_outputs.ingest_outputs_tsv
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

task format_imputation_outputs {
    input {
        String          aggregated_imputation_metrics
        String          chunks_info
        String          failed_chunks
        String          imputed_multisample_vcf
        String          imputed_multisample_vcf_index
        Array[String]?  imputed_single_sample_vcfs
        Array[String]?  imputed_single_sample_vcf_indices
        String          n_failed_chunks
    }

    command <<<
        open_bracket='["'
        close_bracket='"]'

        echo -e "aggregated_imputation_metrics\tchunks_info\tfailed_chunks\tn_failed_chunks\t\
        imputed_multisample_vcf\timputed_multisample_vcf_index\t\
        imputed_single_sample_vcf\timputed_single_sample_vcf_index" \
        > ingestDataset_imputation_outputs.tsv

        echo -e "~{aggregated_imputation_metrics}\t~{chunks_info}\t~{failed_chunks}\t~{n_failed_chunks}\t\
        ~{imputed_multisample_vcf}\t~{imputed_multisample_vcf_index}\t\
        ${open_bracket}~{sep='", "' imputed_single_sample_vcfs}${close_bracket}\t\
        ${open_bracket}~{sep='", "' imputed_single_sample_vcf_indices}${close_bracket}" \
        >> ingestDataset_imputation_outputs.tsv


        python3 << CODE
        import pandas as pd

        tsv_df = pd.read_csv("ingestDataset_imputation_outputs.tsv", sep="\t")
        tsv_df = tsv_df.dropna(axis=1, how="all")  # drop columns if no value (optional outputs etc)

        outputs = tsv_df.to_json("ingestDataset_imputation_outputs.json", orient="records")  # write json file

        CODE
    >>>

    runtime {
        docker: "broadinstitute/horsefish:emerge_scripts"
    }

    output {
        File ingest_outputs_tsv = "ingestDataset_imputation_outputs.tsv"
    }
}

task ingest_outputs_to_tdr {
    input {
        String workspace_name
        String workspace_bucket
        String tdr_dataset_id
        String tdr_target_table_name

        File   outputs_tsv
    }

    command {

        python3 /scripts/emerge/WDL_write_arrays_wdl_outputs_to_TDR_ArraysOutputsTable.py -w ~{workspace_name} \
                                                                          -b ~{workspace_bucket} \
                                                                          -d ~{tdr_dataset_id} \
                                                                          -t ~{tdr_target_table_name} \
                                                                          -f ~{outputs_tsv}
    }

    runtime {
        docker: "broadinstitute/horsefish:emerge_scripts"
    }

    output {
        File ingest_logs = stdout()
    }
}