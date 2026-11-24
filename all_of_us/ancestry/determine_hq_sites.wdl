version 1.0
# Determine which sites in the training data are high quality (hq).
#
# This workflow breaks the detection into two phases.  The first is to process the input VCF, using GATK, and perform
#    all steps, except LD pruning.  Next, this workflow does the LD pruning on the VCF, using Hail.
# This workflow should not need to be run often, as it produces a set of sites to use for training,
#    testing, and applying of ancestry.
# This workflow assumes hg38
#    "High quality sites" are those fulfill the following:
#     AF > 0.1%
#     Biallelic SNPs
#     PASS filtering
#     Have a call rate > 99%
#     LD Pruned
#
# Assumption:  The input vcfs are sharded by chromosome and specified in sorted order.
#           This allows us to run ldPrune in parallel (You're welcome).
#
# Known Issues:
#   This script uses a separate task to do indexing of VCFs.  That is wasteful.  The VCFs could be indexed in the same
#     task that generates them.
#
workflow determine_hq_sites {

    input {
        # VCFs should appear in the order to be merged.  These are training VCFs (e.g. HGDP+1kg)
        Array[File] vcfs

        # A list of info annotations to drop from the training data.  This is useful to manage the size of the VCFs.
        File drop_info_annotations_param_file

        # An identifier for the output set
        String output_prefix

        # Intervals to consider for high quality sites.  This could be a calling interval.
        File intervals

        # Further intersect sites.  Eg, exome intervals.  This can be used to further reduce the area to search for SNPs
        File? intersecting_intervals
    }

    String pipeline_version = "aou_9.0.0"
    String full_vcf_output_prefix = output_prefix + "_full"

    scatter (vcf in vcfs) {
        call simpleFilterForHQSites {
            input:
                vcf_bgz_gts = vcf,
                drop_info_annotations_param_file = drop_info_annotations_param_file,
                intervals = intervals,
                intersecting_intervals = intersecting_intervals
        }

        call ldPrune {
            input:
                vcf_bgz_gts = simpleFilterForHQSites.vcf_bgz_partially_filtered
        }

        call IndexFeatureFile as sitesOnlyIndex {
            input:
                vcf = ldPrune.vcf_sites_only
        }

        call IndexFeatureFileBgz as fullIndex {
            input:
                vcf = ldPrune.vcf_bgz_full
        }
    }

    call MergeVCFs as mergeSitesOnly {
        input:
            input_vcfs = ldPrune.vcf_sites_only,
            input_vcf_indices = sitesOnlyIndex.vcf_idx,
            output_name = output_prefix
    }

    call MergeVCFBgzs as mergeFull {
        input:
            input_vcfs = ldPrune.vcf_bgz_full,
            input_vcf_indices = fullIndex.vcf_idx,
            output_name = full_vcf_output_prefix
    }

    call CreateIntervalListAndUcscBed {
        input :
            vcf=mergeSitesOnly.merged_vcf,
            vcf_idx=mergeSitesOnly.merged_vcf_idx
    }


    output {
        File vcf_bgz_merged_sites_only = mergeSitesOnly.merged_vcf
        File vcf_bgz_merged_sites_only_index = mergeSitesOnly.merged_vcf_idx
        File vcf_bgz_merged_full = mergeFull.merged_vcf
        File vcf_bgz_merged_full_index = mergeFull.merged_vcf_idx
        Array[File] vcf_bgz_fulls = ldPrune.vcf_bgz_full
        Array[File] vcf_bgz_fulls_index = fullIndex.vcf_idx
        File hq_site_ucsc_bed = CreateIntervalListAndUcscBed.ucsc_bed
        File hq_site_interval_list = CreateIntervalListAndUcscBed.interval_list
    }

}
# This task does all of the HQ filtering, in GATK, except LD Pruning.
task simpleFilterForHQSites {

    input {
        File vcf_bgz_gts
        File drop_info_annotations_param_file
        File intervals

        File? intersecting_intervals
    }
    String output_filename = basename(vcf_bgz_gts) + ".partially_filtered.vcf.bgz"

    parameter_meta {
        vcf_bgz_gts: {localization_optional: true}
        intervals: {localization_optional: true}
    }

    command <<<
        set -e
        gatk --java-options "-Xmx4096m" SelectVariants -V ~{vcf_bgz_gts} \
        -L ~{intervals} \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        --exclude-filtered \
        --max-nocall-fraction 0.01 \
        -select "AF>0.001" \
        ~{"-L " + intersecting_intervals + " --interval-set-rule INTERSECTION"} \
        --arguments_file ~{drop_info_annotations_param_file} \
        -O ~{output_filename}
     >>>

    output {
        File vcf_bgz_partially_filtered="~{output_filename}"
    }

    runtime {
        docker:"us.gcr.io/broad-gatk/gatk:4.2.0.0"
        memory: "7 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
    }
}

# Create an LD-pruned VCF file (both full version with genotypes and a sites-only VCF).
# Uses Hail on a single (large) machine.  This is a risk to scaling in the future.
#
# IMPORTANT risks to scaling:
#  - RAM usage numbers are hardcoded.
#  - Hail is being run on a single machine, instead of a Spark cluster.  LD pruning does not
#     parallelize well (w/in a chromosome), but this is still a possible risk.
#
task ldPrune {
    input {
        File vcf_bgz_gts
    }
    String output_filename = basename(vcf_bgz_gts) + ".sites_only.vcf"
    String output_full_filename = basename(vcf_bgz_gts) + ".full.vcf.bgz"

    parameter_meta {
        vcf_bgz_gts: {localization_optional: true}
    }

    command <<<
        set -e
        python3 <<EOF

        import pandas as pd
        import numpy as np
        import hail as hl
        import pyspark

        spark_conf_more_ram = dict()
        spark_conf_more_ram["spark.executor.memory"] = "24g"
        spark_conf_more_ram["spark.driver.memory"] = "24g"
        hl.init(default_reference='GRCh38', idempotent=True, spark_conf=spark_conf_more_ram)

        v = hl.import_vcf("~{vcf_bgz_gts}", force_bgz=True)

        # Perform the ld pruning if there are more than one variants.  If there are 0-1 variants then do not run ldPrune
        # LD Prune may crash if given no variants
        if (v.rows().count() > 1):
            pruned_variant_table = hl.ld_prune(v.GT, r2=0.1, bp_window_size=500000)

            # filter the original MatrixTable down to un-pruned variants
            filtered_ds = v.filter_rows(hl.is_defined(pruned_variant_table[v.row_key]))
        else:
            filtered_ds = v

        # Write a sites-only VCF as output
        hl.export_vcf(filtered_ds.rows(), "~{output_filename}")

        # Write the full VCF
        hl.export_vcf(filtered_ds, "~{output_full_filename}")

        EOF
    >>>

    output {
        File vcf_sites_only="~{output_filename}"
        File vcf_bgz_full="~{output_full_filename}"
    }

    runtime {
        docker: "hailgenetics/hail:0.2.97"
        memory: "123 GB"
        cpu: "4"
        disks: "local-disk 500 HDD"
    }
}

task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_name
    }

    String output_vcf = basename(output_name) + ".vcf.gz"
    String output_vcf_idx = basename(output_vcf) + ".tbi"

    command <<<
        set -e
        gatk --java-options "-Xmx2048m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 700 HDD"
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task MergeVCFBgzs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_name
    }

    String output_vcf = basename(output_name) + ".vcf.bgz"
    String output_vcf_idx = basename(output_vcf) + ".tbi"

    command <<<
                bcftools concat --naive -o ~{output_vcf} ~{sep=" " input_vcfs}
                bcftools index -t ~{output_vcf}
    >>>

    runtime {
        docker: "mgibio/bcftools-cwl:1.12"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 700 HDD"
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task IndexFeatureFile {
    input{
        File vcf
    }
    String output_vcf_idx = basename(vcf) + ".idx"

    command <<<
    set -e
    gatk IndexFeatureFile -I ~{vcf} -O ~{output_vcf_idx}
    >>>

    output {
        File vcf_idx = "~{output_vcf_idx}"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 500 HDD"
    }
}
# TODO: This is a lot of duplicate code just to change the extension on the index file
task IndexFeatureFileBgz {
    input{
        File vcf
    }
    String output_vcf_idx = basename(vcf) + ".tbi"

    command <<<
        set -e
        gatk IndexFeatureFile -I ~{vcf} -O ~{output_vcf_idx}
    >>>

    output {
        File vcf_idx = "~{output_vcf_idx}"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
}

task CreateIntervalListAndUcscBed {
    input {
        File vcf
        File vcf_idx
    }
    String output_vcf = basename(vcf)
    command <<<
        set -e
        gatk VcfToIntervalList -I ~{vcf} -O ~{output_vcf}.interval_list

        gatk IntervalListToBed -I ~{output_vcf}.interval_list -O ~{output_vcf}.ucsc.bed
    >>>

    output {
        File interval_list = "~{output_vcf}.interval_list"
        File ucsc_bed = "~{output_vcf}.ucsc.bed"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.2.6.1"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
}