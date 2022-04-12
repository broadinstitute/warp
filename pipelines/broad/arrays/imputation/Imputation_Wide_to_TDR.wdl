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
        String          output_callset_name
        File            ref_dict
        String          reference_panel_path

        # optional inputs to Imputation.wdl but required for eMerge
        Array[File]     single_sample_vcfs
        Array[File]     single_sample_vcf_indices
        Array[String]   chip_well_barcodes
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

    # call format_imputation_outputs {
    #     input:
    #         imputed_single_sample_vcfs          = Imputation.imputed_single_sample_vcfs,
    #         imputed_single_sample_vcf_indices   = Imputation.imputed_single_sample_vcf_indices,
    #         chip_well_barcodes                  = chip_well_barcodes,
    #         imputed_multisample_vcf             = Imputation.imputed_multisample_vcf,
    #         imputed_multisample_vcf_index       = Imputation.imputed_multisample_vcf_index,
    #         aggregated_imputation_metrics       = Imputation.aggregated_imputation_metrics,
    #         chunks_info                         = Imputation.chunks_info,
    #         failed_chunks                       = Imputation.failed_chunks,
    #         n_failed_chunks                     = Imputation.n_failed_chunks
    # }

    call format_imputation_wide_outputs {
        input:
            imputed_single_sample_vcfs          = Imputation.imputed_single_sample_vcfs,
            imputed_single_sample_vcf_indices   = Imputation.imputed_single_sample_vcf_indices
    }

    # call ingest_outputs_to_tdr {
    #     input:
    #         workspace_name          = workspace_name,
    #         workspace_bucket        = workspace_bucket,
    #         tdr_dataset_id          = tdr_dataset_id,
    #         tdr_target_table_name   = tdr_target_table_name,
    #         outputs_tsv             = format_imputation_outputs.ingest_outputs_tsv
    # }

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
        Array[String]   chip_well_barcodes
        String          n_failed_chunks
    }

    command <<<

        # write header to file
        echo -e "aggregated_imputation_metrics\tchunks_info\tfailed_chunks\tn_failed_chunks\t\
        imputed_multisample_vcf\timputed_multisample_vcf_index\t\
        imputed_single_sample_vcfs\timputed_single_sample_vcf_indices\t\
        chip_well_barcodes" \
        > ingestDataset_imputation_outputs.tsv

        # handle array[type] variables to print as list with double quotes
        imputed_single_sample_vcfs='~{sep='","' imputed_single_sample_vcfs}'
        echo "imputed_single_sample_vcfs"
        echo "[\"${imputed_single_sample_vcfs}\"]"

        imputed_single_sample_vcf_indices='~{sep='","' imputed_single_sample_vcf_indices}'
        echo "imputed_single_sample_vcf_indices"
        echo "[\"${imputed_single_sample_vcf_indices}\"]"

        chip_well_barcodes='~{sep='","' chip_well_barcodes}'
        echo "chip_well_barcodes"
        echo "[\"${chip_well_barcodes}\"]"

        # write file paths to row in tsv file
        echo -e "~{aggregated_imputation_metrics}\t~{chunks_info}\t~{failed_chunks}\t~{n_failed_chunks}\t\
        ~{imputed_multisample_vcf}\t~{imputed_multisample_vcf_index}\t\
        [\"${imputed_single_sample_vcfs}\"]\t\
        [\"${imputed_single_sample_vcf_indices}\"]\t\
        [\"${chip_well_barcodes}\"]" \
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

task format_imputation_wide_outputs{
    input {
        Array[String]?  imputed_single_sample_vcfs
        Array[String]?  imputed_single_sample_vcf_indices
    }

    String prefix="[\""
    String postfix="\"]"

    command <<<

        # handle array[type] variables to print as list with double quotes
        imputed_single_sample_vcfs='~{sep='","' imputed_single_sample_vcfs}'
        imputed_single_sample_vcf_indices='~{sep='","' imputed_single_sample_vcf_indices}'
        export imputed_single_sample_vcfs
        export imputed_single_sample_vcf_indices

        python3 << CODE
        import pandas as pd
        import os

        print("imputed_single_sample_vcfs and indices with os.environ export method")
        os_imputed_single_sample_vcfs = '["' + os.environ["imputed_single_sample_vcfs"] + '"]'
        print("type of imputed_single_sample_vcfs:")
        print(type(os_imputed_single_sample_vcfs))
        
        os_imputed_single_sample_vcf_indices = '["' + os.environ["imputed_single_sample_vcf_indices"] + '"]'
        print("type of imputed_single_sample_vcf_indices:")
        print(type(os_imputed_single_sample_vcf_indices))

        print("imputed_single_sample_vcfs with prefix and postfix method")
        ppt_simple_sample_vcfs=~{prefix}~{sep="\", \"" imputed_single_sample_vcfs}~{postfix}
        print(type(ppt_simple_sample_vcfs))
        print(ppt_simple_sample_vcfs)

        print("imputed_single_sample_vcf_indices with prefix and postfix method")
        ppt_simple_sample_vcf_indices=~{prefix}~{sep="\", \"" imputed_single_sample_vcf_indices}~{postfix}
        print(type(ppt_simple_sample_vcf_indices))
        print(ppt_simple_sample_vcf_indices)

        # tsv_df = pd.DataFrame(columns = ["chip_well_barcode", "imputed_single_sample_vcf", "imputed_single_sample_vcf_index"], sep="\t")
        # sample_dict = {}
        # # for each file in list of imputed vcfs, get chip_well_barcode value
        # for vcf in imputed_single_sample_vcfs:
        #     imputed_vcf_filename = vcf.split("/")[-1]
        #     imputed_vcf_index_filename = imputed_vcf_filename + ".tbi"
        #     chip_well_barcode = filename.split(".")[0]
        #     imputed_vcf_path = vcf
        #     imputed_vcf_index_path = [s for s in imputed_single_sample_vcf_indices if imputed_vcf_index_filename in s][0]

        #     sample_dict["chip_well_barcode"] = chip_well_barcode
        #     sample_dict["imputed_single_sample_vcf"] = imputed_vcf_path
        #     sample_dict["imputed_single_sample_vcf_index"] = imputed_vcf_index_path

        #     tsv_df = tsv_df.append(sample_dict, ignore_index = True)

        # # tsv_df = pd.read_csv("ingestDataset_imputation_wide_outputs.tsv", sep="\t")
        # tsv_df = tsv_df.dropna(axis=1, how="all")  # drop columns if no value (optional outputs etc)
        # outputs = tsv_df.to_json("ingestDataset_imputation_wide_outputs.json", orient="records")  # write json file

        CODE
    >>>

    runtime {
        docker: "broadinstitute/horsefish:emerge_scripts"
    }

    output {
        File ingest_outputs_wide_tsv = "ingestDataset_imputation_wide_outputs.tsv"
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