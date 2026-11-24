version 1.0

# Read in the PCA features from the ancestry table
workflow run_sample_outlier_qc {
    input {
        # This is a tsv of summary stats from GVS.  This is acquired from the GVS team.
        File callset_summary_csv

        # Ancestry results (from run_ancestry.wdl)
        File ancestry_results_tsv

        String output_prefix

        String? metrics_to_check_in
    }

    String metrics_to_check = select_first([metrics_to_check_in, "['snp_count', 'ins_del_ratio', 'del_count', 'ins_count', 'snp_het_homvar_ratio', 'indel_het_homvar_ratio', 'ti_tv_ratio', 'singleton']"])
    String pipeline_version = "aou_9.0.0"

    call join_ancestry_to_stats {
        input:
            ancestry_results_tsv=ancestry_results_tsv,
            callset_summary_csv=callset_summary_csv,
            metrics_to_check=metrics_to_check,
            output_prefix=output_prefix
    }

    call determine_outlier_qc {
        input:
            full_ancestry_ht_tar_gz = join_ancestry_to_stats.full_ancestry_ht_tar_gz,
            metrics_to_check=metrics_to_check,
            input_prefix=output_prefix,
            output_prefix=output_prefix
    }

    output {
        File flagged_samples_tsv = determine_outlier_qc.flagged_samples_tsv
        File all_samples_tsv = determine_outlier_qc.all_samples_tsv
        File ancestry_with_flagged_samples_tar_gz = determine_outlier_qc.all_samples_ht_tar_gz
    }
}

task join_ancestry_to_stats {
    input {
        File ancestry_results_tsv
        File callset_summary_csv
        String metrics_to_check
        String output_prefix
    }
    parameter_meta {
        ancestry_results_tsv: {localization_optional: true}
        callset_summary_csv: {localization_optional: true}
    }
    command <<<
        set -e
        python3 <<EOF

        # Start Hail
        import pandas as pd
        import numpy as np
        import hail as hl
        from bokeh.io import show
        from bokeh.layouts import gridplot
        hl.init(default_reference='GRCh38', idempotent=True)

        # These need to match what is in the callset_summary_tsv
        metrics_to_check = ~{metrics_to_check}

        # Import the tables
        # Note that the keys (research_id and sample_name) of the tables are both being coerced into str
        ancestry_ht = hl.import_table("~{ancestry_results_tsv}", impute=True, key="research_id",
            types={"pca_features":hl.dtype('array<float64>'),
            "probabilities":hl.dtype('array<float64>'),
            "research_id":hl.dtype('str')})
        ancestry_ht = ancestry_ht.rename({"research_id":"s"})

        callset_summary_ht = hl.import_table("~{callset_summary_csv}", impute=True, key="sample_name", delimiter=',',
            types={"sample_name":hl.dtype('str')})

        # Join the metrics from the callset summary to the ancestry table
        for metric in metrics_to_check:
            metric_map = {metric: callset_summary_ht[ancestry_ht.s][metric]}
            ancestry_ht = ancestry_ht.annotate(**metric_map)

        # Sort the sample IDs
        ancestry_ht = ancestry_ht.order_by("s").key_by("s")

        ancestry_ht.write('~{output_prefix}.full_ancestry.ht')

        EOF

        # tar the Hail Tables for output in cromwell
        tar zcvf ~{output_prefix}.full_ancestry.ht.tar.gz ~{output_prefix}.full_ancestry.ht
    >>>

    output {
        File full_ancestry_ht_tar_gz = "~{output_prefix}.full_ancestry.ht.tar.gz"
    }

    runtime {
        docker: "hailgenetics/hail:0.2.67"
        memory: "15 GB"
        cpu: "2"
        disks: "local-disk 500 HDD"
    }
}

task determine_outlier_qc {
    input {
        File full_ancestry_ht_tar_gz
        String metrics_to_check
        String input_prefix
        String output_prefix
    }

    String output_full_tsv_filename = output_prefix + ".full.tsv"
    String output_full_ht_filename = output_prefix + ".full.ht"
    String output_full_ht_tar_gz_filename = output_full_ht_filename + ".tar.gz"
    String output_failed_tsv_filename = output_prefix + ".flagged_samples.tsv"
    String output_failed_ht_filename = output_prefix + ".flagged_samples.ht"
    String output_failed_ht_tar_gz_filename = output_failed_ht_filename + ".tar.gz"

    command <<<
        set -e
        tar zxvf ~{full_ancestry_ht_tar_gz}

        pip install gnomad==0.5.0

        python3 <<EOF
        import pandas as pd
        import numpy as np
        import hail as hl
        import math
        from collections import defaultdict
        from gnomad.sample_qc.filtering import compute_stratified_metrics_filter,compute_qc_metrics_residuals
        hl.init(default_reference='GRCh38', idempotent=True)

        full_ancestry_ht = hl.methods.read_table('~{input_prefix}.full_ancestry.ht')

        metrics_to_check = ~{metrics_to_check}

        qc_metrics = {k:full_ancestry_ht[k] for k in metrics_to_check}

        # Annotate the residuals w/ ancestry to allow stratification
        residuals_ht = compute_qc_metrics_residuals(ht=full_ancestry_ht, pc_scores=full_ancestry_ht.pca_features, qc_metrics=qc_metrics)
        residuals_ht = residuals_ht.annotate(ancestry_pred=full_ancestry_ht[residuals_ht.s].ancestry_pred)
        residuals_ht = residuals_ht.annotate(ancestry_pred_other=full_ancestry_ht[residuals_ht.s].ancestry_pred_other)

        # Prepare data structures for call into gnomad
        #  Default the residual to MADs above and below.  Then override where is makes sense, particularly ones
        #   where we do not want a lower limit.
        metric_threshold_dict = {k + "_residual":(8.0, 8.0) for k in metrics_to_check}
        metric_threshold_dict["singleton_residual"] = (math.inf, 8.0)
        metric_threshold_dict["snp_het_homvar_ratio_residual"] = (math.inf, 8.0)

        qc_metrics_dict = dict(residuals_ht.row_value)
        qc_metrics_dict.pop('ancestry_pred', None)
        qc_metrics_dict.pop('ancestry_pred_other', None)

        # Get the stratified stats for each metric
        strat_ht = compute_stratified_metrics_filter(ht=residuals_ht, qc_metrics=qc_metrics_dict,
            metric_threshold=metric_threshold_dict,
            strata={"pred":residuals_ht.ancestry_pred_other}
        )

        ancestry_strat_ht = full_ancestry_ht.join(strat_ht, how="inner")
        ancestry_strat_ht.write("~{output_full_ht_filename}")
        ancestry_strat_ht.export("~{output_full_tsv_filename}")

        flagged_sample_list_ht = ancestry_strat_ht.filter(ancestry_strat_ht.qc_metrics_filters.size() != 0)
        flagged_sample_list_ht.write("~{output_failed_ht_filename}")
        flagged_sample_list_ht.export("~{output_failed_tsv_filename}")
        EOF

        # tar the Hail Tables for output in cromwell
        tar zcvf ~{output_full_ht_tar_gz_filename} ~{output_full_ht_filename}
        tar zcvf ~{output_failed_ht_tar_gz_filename} ~{output_failed_ht_filename}
    >>>

    output {
        File flagged_samples_tsv = "~{output_failed_tsv_filename}"
        File flagged_samples_ht_tar_gz = "~{output_failed_ht_tar_gz_filename}"
        File all_samples_tsv = "~{output_full_tsv_filename}"
        File all_samples_ht_tar_gz = "~{output_full_ht_tar_gz_filename}"
    }

    runtime {
        docker: "hailgenetics/hail:0.2.67"
        memory: "15 GB"
        cpu: "4"
        disks: "local-disk 500 HDD"
    }
}



