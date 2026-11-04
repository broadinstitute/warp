version 1.0

workflow run_sample_outlier_qc_plotting {
    input {
        File aou_demographics_tsv

        String output_prefix

        # Must have been created from run_sample_outlier_qc.wdl.  Input prefix was the output prefix specified for that
        #  workflow
        File ancestry_with_flagged_samples_tar_gz
        String input_prefix
    }
    String pipeline_version="aou_9.0.0"

    call join_ancestry_to_demographics {
        input:
            aou_demographics_tsv=aou_demographics_tsv,
            output_prefix=output_prefix,
            ancestry_with_flagged_samples_tar_gz=ancestry_with_flagged_samples_tar_gz,
            input_prefix=input_prefix
    }

    call plot_first_pcs {
        input:
            ancestry_with_flagged_samples_demographics_tar_gz=join_ancestry_to_demographics.ancestry_with_flagged_samples_demographics_tar_gz,
            aou_demographics_tsv=aou_demographics_tsv,
            input_prefix=input_prefix,
            output_prefix=output_prefix
    }

    call plot_metrics_and_fitting {
        input:
            ancestry_with_flagged_samples_demographics_tar_gz=join_ancestry_to_demographics.ancestry_with_flagged_samples_demographics_tar_gz,
            aou_demographics_tsv=aou_demographics_tsv,
            input_prefix=input_prefix,
            output_prefix=output_prefix
    }

    output {
        File pcs_plot = plot_first_pcs.pcs_plot
        File metrics_plot = plot_metrics_and_fitting.metrics_plot
    }
}

task join_ancestry_to_demographics {
    input {
        File aou_demographics_tsv
        String output_prefix
        File ancestry_with_flagged_samples_tar_gz
        String input_prefix
    }
    command <<<
        set -e

        tar zxvf ~{ancestry_with_flagged_samples_tar_gz}

        python3 <<EOF

        # Start Hail
        import pandas as pd
        import numpy as np
        import hail as hl
        from bokeh.io import show
        from bokeh.layouts import gridplot
        from bokeh.io import output_file
        hl.init(default_reference='GRCh38', idempotent=True)

        # Load the ancestry information
        ancestry_strat_ht = hl.methods.read_table('~{input_prefix}.full.ht')

        #  Convert the ID into a string and set that as the key field.
        ancestry_strat_key_field = "s"
        ancestry_strat_ht = ancestry_strat_ht.annotate(s_str=hl.str(ancestry_strat_ht[ancestry_strat_key_field]))
        ancestry_strat_ht = ancestry_strat_ht.key_by('s_str')
        ancestry_strat_ht = ancestry_strat_ht.transmute(s=hl.str(ancestry_strat_ht[ancestry_strat_key_field]))
        ancestry_strat_ht = ancestry_strat_ht.key_by(ancestry_strat_key_field)
        ancestry_strat_ht.describe()

        # Load demo graphic information.  Assumes the delimieter is '\t' and that the esacpe of cells with special chars
        #  is a double quote.
        # We convert the research ID into a str.  We might be better off changing the key to an int in the future, but
        #  that would require an upstream change.
        demographics_ht = hl.methods.import_table("~{aou_demographics_tsv}", key="research_id", impute=True, quote='"', types={"research_id":hl.tstr})
        demographics_ht = demographics_ht.rename({"research_id":ancestry_strat_key_field})
        demographics_ht.describe()

        # Read in race/ethnicity from demographics file (field Race_WhatRaceEthnicity)
        ancestry_strat_ht = ancestry_strat_ht.annotate(race=demographics_ht[ancestry_strat_ht[ancestry_strat_key_field]]['Race_WhatRaceEthnicity'])
        ancestry_strat_ht.write('~{output_prefix}.ancestry_with_flagged_samples_demographics.ht')

        EOF

        tar zcvf ~{output_prefix}.ancestry_with_flagged_samples_demographics.ht.tar.gz ~{output_prefix}.ancestry_with_flagged_samples_demographics.ht
    >>>

    output {
        File ancestry_with_flagged_samples_demographics_tar_gz = "~{output_prefix}.ancestry_with_flagged_samples_demographics.ht.tar.gz"
    }
    runtime {
        docker: "hailgenetics/hail:0.2.67"
        memory: "7 GB"
        cpu: "4"
        disks: "local-disk 100 HDD"
    }
}

task plot_first_pcs {
    input {
        File ancestry_with_flagged_samples_demographics_tar_gz
        File aou_demographics_tsv
        String input_prefix
        String output_prefix
    }

    command <<<
        set -e

        tar zxvf ~{ancestry_with_flagged_samples_demographics_tar_gz}

        python3 <<EOF

        # Start Hail
        import pandas as pd
        import numpy as np
        import hail as hl
        from bokeh.io import show
        from bokeh.layouts import gridplot
        from bokeh.io import output_file
        hl.init(default_reference='GRCh38', idempotent=True)

        ancestry_strat_ht = hl.methods.read_table('~{input_prefix}.ancestry_with_flagged_samples_demographics.ht')
        output_file("~{output_prefix}.pc1vspc2.html")
        scatter_plot = hl.plot.scatter(ancestry_strat_ht.pca_features[0], ancestry_strat_ht.pca_features[1], label=ancestry_strat_ht['ancestry_pred_other'], hover_fields={'self_reported_race':ancestry_strat_ht['race'].replace('WhatRaceEthnicity_', ''), 'rid': ancestry_strat_ht['s']}, title="~{output_prefix} (PCs 1 and 2)")
        show(scatter_plot)
        EOF

    >>>

    output {
        File pcs_plot = "~{output_prefix}.pc1vspc2.html"
    }

    runtime {
        docker: "hailgenetics/hail:0.2.67"
        memory: "7 GB"
        cpu: "4"
        disks: "local-disk 100 HDD"
    }
}

task plot_metrics_and_fitting {
    input {
        File ancestry_with_flagged_samples_demographics_tar_gz
        File aou_demographics_tsv
        String input_prefix
        String output_prefix
    }

    command <<<
        set -e

        tar zxvf ~{ancestry_with_flagged_samples_demographics_tar_gz}

        python3 <<EOF

        # Start Hail
        import pandas as pd
        import numpy as np
        import hail as hl
        from bokeh.io import show
        from bokeh.layouts import gridplot
        from bokeh.io import output_file
        from bokeh.plotting import figure
        from bokeh.models.widgets import Panel, Tabs

        spark_conf_more_ram = dict()
        spark_conf_more_ram["spark.executor.memory"] = "4g"
        spark_conf_more_ram["spark.executor.cores"] = "4"
        spark_conf_more_ram["spark.driver.memory"] = "4g"
        hl.init(default_reference='GRCh38', idempotent=True, spark_conf=spark_conf_more_ram)

        ancestry_strat_ht = hl.methods.read_table('~{input_prefix}.ancestry_with_flagged_samples_demographics.ht')

        def create_metric_panel(metric:str, ancestry_ht:hl.Table, pca_ind:int = 0):

            tmp = ancestry_ht.pca_features.collect()
            tmp2 = ancestry_ht.lms[metric].beta.collect()[0]
            tmp3 = ancestry_ht.filter(ancestry_strat_ht["fail_" + metric + "_residual"]).pca_features.collect()
            tmp4 = ancestry_ht.filter(ancestry_strat_ht["fail_" + metric + "_residual"])[metric].collect()
            tmp5 = ancestry_ht.filter(ancestry_strat_ht["fail_" + metric + "_residual"]).s.collect()

            xvals = [x[pca_ind] for x in tmp]
            xvals.sort()
            yvals = [tmp2[pca_ind + 17]*(x**2) + x*tmp2[pca_ind + 1] + tmp2[0] for x in xvals]

            xvals_fail = [x[pca_ind] for x in tmp3]
            yvals_fail = [y for y in tmp4]
            labels_fail = [s for s in tmp5]

            # Pretty print metrics
            pretty_metrics = {
                'snp_count':"SNP count",
                'ins_del_ratio': "Insertion-deletion ratio",
                'del_count': "Deletion count",
                'ins_count': "Insertion count",
                'snp_het_homvar_ratio': "SNP het-homvar ratio",
                'indel_het_homvar_ratio': "Indel het-homvar ratio",
                'ti_tv_ratio': "TiTv ratio",
                'singleton': "Variant count not in gnomAD"
            }


            scatter_plot = hl.plot.scatter(ancestry_ht.pca_features[pca_ind], ancestry_ht[metric],
                label=ancestry_ht['ancestry_pred_other'], hover_fields={'self_reported_race':ancestry_ht['race'].replace('WhatRaceEthnicity_', ''), 'rid': ancestry_ht['s']},
                title="(PC " + str(pca_ind+1) + " vs " + metric +") " + str(len(xvals_fail)) + " failures",
                xlabel="PC " + str(pca_ind+1), ylabel=pretty_metrics.get(metric, metric))

            scatter_plot.xaxis.axis_label_text_font_size  = '30px'
            scatter_plot.yaxis.axis_label_text_font_size  = '30px'
            scatter_plot.xaxis.major_label_text_font_size = "15px"
            scatter_plot.yaxis.major_label_text_font_size = "15px"

            scatter_plot.line(xvals, yvals, line_width=2)
            c1 = scatter_plot.circle(xvals_fail, yvals_fail, size=20, fill_color=None, alpha=0.25)
            return scatter_plot

        pca_ind = 0

        # Set up the file for output
        output_file("~{output_prefix}.metrics.html")

        metrics_to_check = ['snp_count', 'ins_del_ratio', 'del_count', 'ins_count', 'snp_het_homvar_ratio', 'indel_het_homvar_ratio', 'ti_tv_ratio', 'singleton']
        panels = [Panel(child=create_metric_panel(m, ancestry_strat_ht, pca_ind), title=m + " and PC " + str(pca_ind+1)) for m in metrics_to_check]
        tabs = Tabs(tabs=panels)
        show(tabs)

        EOF

    >>>

    output {
        File metrics_plot = "~{output_prefix}.metrics.html"
    }

    runtime {
        docker: "hailgenetics/hail:0.2.67"
        memory: "26 GB"
        cpu: "4"
        disks: "local-disk 100 HDD"
    }
}