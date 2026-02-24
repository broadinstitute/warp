version 1.0

import "https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/v2.5_MongoSwirl_Single/MongoTasks_v2_5_Single.wdl" as MongoTasks_Single

workflow AlignAndCallR2 {
    meta {
        description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
    }

    input {
        File unmapped_bam
        String sample_name
        String suffix

        File mt_interval_list

        File mt_self
        File mt_self_index
        File mt_self_dict

        File self_cat
        File self_cat_index
        File self_cat_dict

        File mt_self_shifted
        File mt_self_shifted_index
        File mt_self_shifted_dict

        File self_shifted_cat
        File self_shifted_cat_index
        File self_shifted_cat_dict

        File blacklisted_sites
        File blacklisted_sites_index

        File force_call_vcf
        File force_call_vcf_idx
        File force_call_vcf_shifted
        File force_call_vcf_shifted_idx

        File shift_back_chain

        File non_control_interval
        File control_shifted

        File? gatk_override
        String? gatk_docker_override
        String gatk_version = "4.2.6.0"
        String? m2_extra_args
        String? m2_filter_extra_args
        Float? vaf_filter_threshold
        Float? f_score_beta
        Boolean compress_output_vcf

        Float? verifyBamID
        String hasContamination
        Float contamination_major
        Float contamination_minor

        # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
        # affected by this number. Default is 151.
        Int? max_read_length

        #Optional runtime arguments
        Int? preemptible_tries
        Int? n_cpu
        Int? n_cpu_bwa
    }

    parameter_meta {
        unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
    }

    call MongoTasks_Single.MongoAlignToMtRegShiftedAndMetrics as AlignToMtRegShiftedAndMetrics {
        input:
            input_bam = unmapped_bam,
            sample_base_name = sample_name,
            suffix = suffix,

            mt = mt_self,
            mt_index = mt_self_index,
            mt_dict = mt_self_dict,

            mt_cat = self_cat,
            mt_cat_index = self_cat_index,
            mt_cat_dict = self_cat_dict,

            mt_shifted = mt_self_shifted,
            mt_shifted_index = mt_self_shifted_index,
            mt_shifted_dict = mt_self_shifted_dict,

            mt_shifted_cat = self_shifted_cat,
            mt_shifted_cat_index = self_shifted_cat_index,
            mt_shifted_cat_dict = self_shifted_cat_dict,

            mt_interval_list = mt_interval_list,

            read_length = max_read_length,
            coverage_cap = 100000,

            preemptible_tries = preemptible_tries,
            n_cpu = n_cpu_bwa
    }

    Int M2_mem = if AlignToMtRegShiftedAndMetrics.mean_coverage > 25000 then 14 else 7

    call MongoTasks_Single.MongoCallMtAndShifted as CallMtAndShifted {
        input:
            sample_base_name = sample_name,
            suffix = suffix,
        # Everything is called except the control region.
            input_bam = AlignToMtRegShiftedAndMetrics.mt_aligned_bam,
            input_bai = AlignToMtRegShiftedAndMetrics.mt_aligned_bai,
            mt_self = mt_self,
            mt_self_index = mt_self_index,
            mt_self_dict = mt_self_dict,
            mt_interval_list = non_control_interval,
            m2_extra_args = select_first([m2_extra_args, ""]),
            force_call_vcf = force_call_vcf,
            force_call_vcf_idx = force_call_vcf_idx,

        # Only the control region is now called.
            shifted_input_bam = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bam,
            shifted_input_bai = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bai,
            shifted_mt_self = mt_self_shifted,
            shifted_mt_self_index = mt_self_shifted_index,
            shifted_mt_self_dict = mt_self_shifted_dict,
            shifted_mt_interval_list = control_shifted,
            shifted_m2_extra_args = select_first([m2_extra_args, ""]),
            shifted_force_call_vcf = force_call_vcf_shifted,
            shifted_force_call_vcf_idx = force_call_vcf_shifted_idx,

            compress = compress_output_vcf,
            gatk_override = gatk_override,
            gatk_docker_override = gatk_docker_override,
            gatk_version = gatk_version,
            mem = M2_mem,
            preemptible_tries = preemptible_tries,
            n_cpu = n_cpu
    }

    call MongoTasks_Single.MongoLiftoverCombineMergeFilterContamSplit as LiftoverCombineMergeFilterContamSplit {
        input:
            sample_base_name = sample_name,
            suffix = suffix,

            mt_self = mt_self,
            mt_self_index = mt_self_index,
            mt_self_dict = mt_self_dict,
            shifted_vcf = CallMtAndShifted.shifted_raw_vcf,
            shifted_vcf_idx = CallMtAndShifted.shifted_raw_vcf_idx,
            non_shifted_vcf = CallMtAndShifted.raw_vcf,
            non_shifted_vcf_idx = CallMtAndShifted.raw_vcf_idx,
            shifted_stats = CallMtAndShifted.shifted_stats,
            non_shifted_stats = CallMtAndShifted.stats,
            shift_back_chain = shift_back_chain,
            blacklisted_sites = blacklisted_sites,
            blacklisted_sites_index = blacklisted_sites_index,

            hasContamination = hasContamination,
            contamination_major = contamination_major,
            contamination_minor = contamination_minor,
            verifyBamID = verifyBamID,

            compress = compress_output_vcf,
            m2_extra_filtering_args = m2_filter_extra_args,
            max_alt_allele_count = 4,
            vaf_filter_threshold = vaf_filter_threshold,
            f_score_beta = f_score_beta,

            gatk_override = gatk_override,
            gatk_docker_override = gatk_docker_override,
            gatk_version = gatk_version,
            preemptible_tries = preemptible_tries
    }

    output {
        File mt_aligned_bam = AlignToMtRegShiftedAndMetrics.mt_aligned_bam
        File mt_aligned_bai = AlignToMtRegShiftedAndMetrics.mt_aligned_bai
        File mt_aligned_shifted_bam = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bam
        File mt_aligned_shifted_bai = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bai
        File nuc_mt_aligned_bam = AlignToMtRegShiftedAndMetrics.nuc_and_mt_aligned_bam
        File nuc_mt_aligned_bai = AlignToMtRegShiftedAndMetrics.nuc_and_mt_aligned_bai
        File nuc_mt_shifted_aligned_bam = AlignToMtRegShiftedAndMetrics.nuc_and_shifted_mt_aligned_bam
        File nuc_mt_shifted_aligned_bai = AlignToMtRegShiftedAndMetrics.nuc_and_shifted_mt_aligned_bai
        File out_vcf = LiftoverCombineMergeFilterContamSplit.filtered_vcf
        File out_vcf_idx = LiftoverCombineMergeFilterContamSplit.filtered_vcf_idx
        File split_vcf = LiftoverCombineMergeFilterContamSplit.split_vcf
        File split_vcf_idx = LiftoverCombineMergeFilterContamSplit.split_vcf_idx
        File duplicate_metrics = AlignToMtRegShiftedAndMetrics.duplicate_metrics
        File coverage_metrics = AlignToMtRegShiftedAndMetrics.wgs_metrics
        File theoretical_sensitivity_metrics = AlignToMtRegShiftedAndMetrics.theoretical_sensitivity
        Int mean_coverage = AlignToMtRegShiftedAndMetrics.mean_coverage
        Float median_coverage = AlignToMtRegShiftedAndMetrics.median_coverage
    }
}