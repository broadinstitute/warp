version 1.0

# This workflow will take a VDS then:
#  1. Split it into separate processing for each chromsome.
#  2. filter it by a BED file
#  3. Write out as a VCF (.bgz)
#  4. Write out as a sites-only VCF (.bgz)
#
# No merging of the VCFs occur
# This workflow assumes hg38
# Does not work with requester pays buckets in cromwell
workflow vds_to_vcf {
    input {
        String vds_gs_url

        String bed_gs_url

        # This is the number of partitions to apply to the entire VDS, not each contig.
        # Please note that 70000 is what we use for big VDS, such as the 245k AoU WGS callset
        Int n_partitions = 2000

        String output_prefix

        # This should be ordered
        Array[String] contigs
    }
    String pipeline_version = "aou_9.0.0"

    scatter (contig in contigs) {
        call process_vds {
            input:
                vds_gs_url=vds_gs_url,
                bed_gs_url=bed_gs_url,
                chromosome=contig,
                n_parts=n_partitions,
                output_prefix=output_prefix
        }
    }

    call create_fofn {
        input:
            output_prefix=output_prefix,
            file_urls1=process_vds.vcf,
            file_urls2=process_vds.vcf_idx
    }

    output {
        Array[File] vcfs = process_vds.vcf
        Array[File] vcfs_tbis = process_vds.vcf_idx
        Array[File] vcfs_so = process_vds.vcf_so
        Array[File] vcfs_so_tbis = process_vds.vcf_so_idx
        File vcfs_list = create_fofn.fofn1
        File vcfs_idx_list = create_fofn.fofn2
    }
}

task process_vds {
    input {
        String vds_gs_url
        String bed_gs_url
        String chromosome
        Int n_parts
        String output_prefix
        String docker = "hailgenetics/hail:0.2.134-py3.11"
        Int memory_gb = 624
        Int cpu = 96
        Int bootDiskSize_gb = 500
        Int disk_gb = 1000
    }
    command <<<
        set -e

        # Remove this when you are running on a Spark cluster.
        mkdir -p tmp_dir
        export TMPDIR=${PWD}/tmp_dir/

        python3 <<EOF
        import datetime
        def print_date():
            now = datetime.datetime.now()
            print(str(now))

        # Start Hail
        import pandas as pd
        import numpy as np
        import hail as hl

        spark_conf_more_ram = dict()
        spark_conf_more_ram["spark.executor.memory"] = "12g"
        spark_conf_more_ram["spark.executor.cores"] = "4"
        spark_conf_more_ram["spark.driver.memory"] = "128g"
        hl.init(default_reference='GRCh38', idempotent=False, spark_conf=spark_conf_more_ram)

        print(hl.spark_context().master)

        # Repartition the VDS.
        print("Repartitioning the VDS")
        print_date()
        vds1 = hl.vds.read_vds("~{vds_gs_url}")
        partition_intervals = vds1.reference_data._calculate_new_partitions(~{n_parts})
        callset = hl.vds.read_vds("~{vds_gs_url}", intervals=partition_intervals)
        print_date()

        # Create the interval table for filtering
        print("Creating interval table")
        print_date()
        interval_table = hl.import_bed("~{bed_gs_url}", reference_genome='GRCh38')
        print_date()

        # Filter VDS by chromosome
        callset = hl.vds.filter_chromosomes(callset, keep=["~{chromosome}"])

        # Filter VDS by the bed file.  Note that we will also be updating the reference data to help with processing
        #  downstream
        print("Filtering by intervals...")
        print_date()
        callset = hl.vds.filter_intervals(callset, interval_table, split_reference_blocks=True)
        print_date()

        # Densify
        # Update the variant data to include GT and convert FT to a string.
        #  Also, handle the case where we have no FT.
        vd_gt = callset.variant_data.transmute_entries(
            GT = hl.vds.lgt_to_gt(callset.variant_data.LGT, callset.variant_data.LA)
        )
        if 'FT' in vd_gt.entry:
            vd_gt = vd_gt.transmute_entries(FT = hl.if_else(vd_gt.FT, "PASS", "FAIL"))

        # Drop gvcf_info (if it exists) since it causes issues in test data we have seen.
        if 'gvcf_info' in vd_gt.entry:
            vd_gt = vd_gt.drop('gvcf_info')

        d_callset = hl.vds.to_dense_mt(hl.vds.VariantDataset(callset.reference_data, vd_gt))
        d_callset = d_callset.annotate_rows(gt_stats = hl.agg.call_stats(d_callset.GT, d_callset.alleles))
        d_callset = d_callset.rename({"gt_stats":"info"})
        d_callset.describe()

        # Prune the first entry in AC and AF.  Hail has a convention of including the reference info, which is
        #  against convention when rendered into a VCF
        d_callset = d_callset.transmute_rows(info=d_callset.info.annotate(AF=d_callset.info.AF[1:], AC=d_callset.info.AC[1:]))
        d_callset.describe()

        # Convert to VCF
        print("Doing entire pipeline through full VCFs")
        print_date()
        hl.export_vcf(d_callset, "~{output_prefix}.~{chromosome}.vcf.bgz", tabix=True)
        print_date()
        print("Writing sites-only VCFs, but using a pseudo checkpoint")
        print_date()
        mt = hl.import_vcf("~{output_prefix}.~{chromosome}.vcf.bgz", min_partitions=96)
        hl.export_vcf(mt.rows(), "~{output_prefix}.~{chromosome}.so.vcf.bgz", tabix=True)
        print_date()

        EOF
    >>>
    output {
        File vcf = "~{output_prefix}.~{chromosome}.vcf.bgz"
        File vcf_idx = "~{output_prefix}.~{chromosome}.vcf.bgz.tbi"
        File vcf_so = "~{output_prefix}.~{chromosome}.so.vcf.bgz"
        File vcf_so_idx = "~{output_prefix}.~{chromosome}.so.vcf.bgz.tbi"
    }
    runtime {
        docker: docker
        memory: "~{memory_gb} GiB"
        cpu: cpu
        disks: "local-disk ~{disk_gb} HDD"
        bootDiskSizeGb: bootDiskSize_gb
    }
}

task create_fofn {
    input {
        Array[String] file_urls1
        Array[String] file_urls2
        String output_prefix
        String docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
        Int memory_gb = 3
        Int cpu = 1
        Int disk_gb = 100
    }
    File fofn1_in = write_lines(file_urls1)
    File fofn2_in = write_lines(file_urls2)
    command <<<
        set -e
        cp ~{fofn1_in} ~{output_prefix}.fofn1.txt
        cp ~{fofn2_in} ~{output_prefix}.fofn2.txt
    >>>
    output {
        File fofn1 = "~{output_prefix}.fofn1.txt"
        File fofn2 = "~{output_prefix}.fofn2.txt"
    }
    runtime {
        docker: docker
        memory: "~{memory_gb} GiB"
        cpu: cpu
        disks: "local-disk ~{disk_gb} HDD"
    }
}