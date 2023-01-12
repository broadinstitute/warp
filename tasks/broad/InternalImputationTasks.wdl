version 1.0

task FormatImputationOutputs {
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
        imputed_single_sample_vcf_indices='~{sep='","' imputed_single_sample_vcf_indices}'
        chip_well_barcodes='~{sep='","' chip_well_barcodes}'

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
        docker: "gcr.io/emerge-production/emerge_wdls:v.1.0"
    }

    output {
        File ingest_outputs_tsv = "ingestDataset_imputation_outputs.tsv"
        File ingest_outputs_json = "ingestDataset_imputation_outputs.json"
    }
}

task FormatImputationWideOutputs{
    input {
        Array[String]?  imputed_single_sample_vcfs
        Array[String]?  imputed_single_sample_vcf_indices
    }

    command <<<

        python3 << CODE
        import pandas as pd
        import os

        print("imputed_single_sample_vcfs")
        single_sample_vcfs = [ x for x in [ "~{sep='", "' imputed_single_sample_vcfs}" ]  if x != "" ]

        print("imputed_single_sample_vcf_indices")
        single_sample_vcf_indices = [ x for x in [ "~{sep='", "' imputed_single_sample_vcf_indices}" ]  if x != "" ]

        print("creating dataframe")
        all_samples = []

        print("getting vcf + vcf index file names and paths and chipwell barcode")
        # for each file in list of imputed vcfs, get chip_well_barcode value
        for vcf in single_sample_vcfs:
            sample_dict = {}
            imputed_vcf_filename = vcf.split("/")[-1]
            imputed_vcf_index_filename = imputed_vcf_filename + ".tbi"
            chip_well_barcode = imputed_vcf_filename.split(".")[0]
            imputed_vcf_path = vcf
            imputed_vcf_index_path = [s for s in single_sample_vcf_indices if imputed_vcf_index_filename in s][0]

            sample_dict["chip_well_barcode"] = chip_well_barcode
            sample_dict["imputed_single_sample_vcf"] = imputed_vcf_path
            sample_dict["imputed_single_sample_vcf_index"] = imputed_vcf_index_path

            print("single vcf dictionary for vcf with name:" + imputed_vcf_filename)
            print(sample_dict)

            print("appending single vcf dict to dataframe")
            all_samples.append(sample_dict)
            print("list of dictionaries after adding" + chip_well_barcode)
            print(all_samples)

        print("writing final dataframe to tsv and json file")
        tsv_df = pd.DataFrame(all_samples)
        print("dataframe after adding all imputation samples")
        tsv_df = tsv_df.dropna(axis=1, how="all")  # drop columns if no value (optional outputs etc)

        # write dataframe to tsv
        tsv_df.to_csv("ingestDataset_imputation_wide_outputs.tsv", index=False , sep="\t")
        print("finished writing dataframe of split out imputation outputs to tsv file")
        # write dataframe to json
        outputs = tsv_df.to_json("ingestDataset_imputation_wide_outputs.json", orient="records")  # write json file

        CODE
    >>>

    runtime {
        docker: "gcr.io/emerge-production/emerge_wdls:v.1.0"
    }

    output {
        File ingest_outputs_wide_tsv = "ingestDataset_imputation_wide_outputs.tsv"
        File ingest_outputs_wide_json = "ingestDataset_imputation_wide_outputs.json"
    }
}

task TriggerPrsWithImputationTsv {
    input {
        File    run_task
        File    imputation_outputs_tsv
        String  trigger_bucket_path
        String  timestamp
    }

    command {
        destination_file_name=~{timestamp}"_ingestDataset_imputation_outputs.tsv"
        gsutil cp ~{imputation_outputs_tsv} ~{trigger_bucket_path}$destination_file_name
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    }

    output {
        File trigger_prs_cf_log = stdout()
    }
}