version 1.0
# BroadInternalArrays is now deprecated 2025-03-06
import "../../../../../pipelines/wdl/arrays/single_sample/Arrays.wdl" as ArraysPipeline
import "../../../../../tasks/wdl/InternalArraysTasks.wdl" as InternalArraysTasks
import "../../../../../tasks/wdl/InternalTasks.wdl" as InternalTasks

workflow BroadInternalArrays {
    meta {
        description: "Push outputs of Arrays.wdl to TDR dataset table ArraysOutputsTable."
    }

    String pipeline_version = "1.1.15"

    input {
        # inputs to wrapper task
        String  tdr_dataset_id
        String  tdr_target_table_name

        # required inputs to Arrays.wdl
        String  chip_well_barcode
        String  sample_alias
        String  sample_lsid
        String  reported_gender
        String  lab_batch
        File    red_idat_cloud_path
        File    green_idat_cloud_path
        File    ref_fasta
        File    ref_fasta_index
        File    ref_dict
        File    dbSNP_vcf
        File    dbSNP_vcf_index
        File    haplotype_database_file
        File    variant_rsids_file
        Int     disk_size
        Int     preemptible_tries
        String  environment
        File    vault_token_path
    }

    call ArraysPipeline.Arrays {
        input:
            chip_well_barcode          = chip_well_barcode,
            sample_alias               = sample_alias,
            sample_lsid                = sample_lsid,
            reported_gender            = reported_gender,
            red_idat_cloud_path        = red_idat_cloud_path,
            green_idat_cloud_path      = green_idat_cloud_path, 
            ref_fasta                  = ref_fasta, 
            ref_fasta_index            = ref_fasta_index, 
            ref_dict                   = ref_dict, 
            dbSNP_vcf                  = dbSNP_vcf, 
            dbSNP_vcf_index            = dbSNP_vcf_index, 
            haplotype_database_file    = haplotype_database_file,
            variant_rsids_file         = variant_rsids_file, 
            disk_size                  = disk_size, 
            preemptible_tries          = preemptible_tries, 
            environment                = environment, 
            vault_token_path           = vault_token_path
    }

    call InternalArraysTasks.FormatArraysOutputs {
        input:
            chip_well_barcode_output                            = Arrays.chip_well_barcode_output,
            analysis_version_number_output                      = Arrays.analysis_version_number_output,
            gtc_file                                            = Arrays.gtc_file,
            output_vcf                                          = Arrays.output_vcf,
            output_vcf_index                                    = Arrays.output_vcf_index,
            baf_regress_metrics_file                            = Arrays.baf_regress_metrics_file,
            arrays_variant_calling_detail_metrics_file          = Arrays.arrays_variant_calling_detail_metrics_file,
            arrays_variant_calling_summary_metrics_file         = Arrays.arrays_variant_calling_summary_metrics_file,
            arrays_variant_calling_control_metrics_file         = Arrays.arrays_variant_calling_control_metrics_file,
            fingerprint_detail_metrics_file                     = Arrays.fingerprint_detail_metrics_file,
            fingerprint_summary_metrics_file                    = Arrays.fingerprint_summary_metrics_file,
            genotype_concordance_summary_metrics_file           = Arrays.genotype_concordance_summary_metrics_file,
            genotype_concordance_detail_metrics_file            = Arrays.genotype_concordance_detail_metrics_file,
            genotype_concordance_contingency_metrics_file       = Arrays.genotype_concordance_contingency_metrics_file,
            lab_batch                                           = lab_batch
    }

    call InternalTasks.IngestOutputsToTDR {
        input:
            tdr_dataset_id          = tdr_dataset_id,
            tdr_target_table_name   = tdr_target_table_name,
            outputs_tsv             = FormatArraysOutputs.ingest_outputs_tsv
    }

    output {
        String chip_well_barcode_output = Arrays.chip_well_barcode_output
        Int analysis_version_number_output = Arrays.analysis_version_number_output
        File gtc_file = Arrays.gtc_file
        File? output_vcf = Arrays.output_vcf
        File? output_vcf_index = Arrays.output_vcf_index
        File? baf_regress_metrics_file = Arrays.baf_regress_metrics_file
        File arrays_variant_calling_detail_metrics_file = Arrays.arrays_variant_calling_detail_metrics_file
        File? arrays_variant_calling_summary_metrics_file = Arrays.arrays_variant_calling_summary_metrics_file
        File? arrays_variant_calling_control_metrics_file = Arrays.arrays_variant_calling_control_metrics_file
        File? fingerprint_detail_metrics_file = Arrays.fingerprint_detail_metrics_file
        File? fingerprint_summary_metrics_file = Arrays.fingerprint_summary_metrics_file
        File? genotype_concordance_summary_metrics_file = Arrays.genotype_concordance_summary_metrics_file
        File? genotype_concordance_detail_metrics_file  = Arrays.genotype_concordance_detail_metrics_file
        File? genotype_concordance_contingency_metrics_file = Arrays.genotype_concordance_contingency_metrics_file
    }
}