## SET UP THE ENVIRONMENT
from datetime import datetime
import os
import pandas as pd
import numpy as np
import argparse
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--contig',                     help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',                help='GCS bucket used for spark.local.dir, likely the cluster temp bucket')
# INPUTS
parser.add_argument('--input_aou_vds_url',          help='GCS url to the AoU-only VDS file to process (e.g. gs://bucket/path/to/file.vds)')
# OUTPUTS
parser.add_argument('--output_aou_vcf_url',         help='GCS url where the AoU-only VCF output will be copied')
parser.add_argument('--output_aou_vcf_header_url',  help='GCS url where the AoU-only VCF header output will be copied')
parser.add_argument('--output_report_url',          help='GCS url where the report output will be copied')
args = parser.parse_args()


def run_subprocess(cmd, errorMessage):
    if isinstance(cmd, list):
        cmd = ' '.join(cmd)
    try:
        # print("running command: " + cmd)
        return subprocess.check_output(
            cmd, shell=True, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        print(errorMessage)
        print("Exited with " + str(e.returncode) + "-" + e.output)
        exit(1)


## START THE TIME COUNTER
start = datetime.now()

## DEFINE BUCKETS
tmp_bucket = args.temp_bucket

## RIG DETAILS AND COST
# Main:
# CPUs = 32
# Memory = 208
# Disk size = 150
# Workers:
# Workers = 2
# Preemptibles = 50
# CPUs = 4
# Memory = 26
# Disk size = 150
rig_cost=7.80

## SET UP HAIL ENVIRONMENT (+ PLOTTING)
import hail as hl

## BOOT UP HAIL
temp_dir = f'{tmp_bucket}/Stage_1/temp'

spark_conf_more_ram = dict()
spark_conf_more_ram["spark.executor.memory"] = "8g"
spark_conf_more_ram["spark.executor.cores"] = "4"
spark_conf_more_ram["spark.driver.memory"] = "60g"
spark_conf_more_ram["spark.driver.cores"] = "32"
spark_conf_more_ram["spark.local.dir"] = temp_dir

hl.init(idempotent=True,
        tmp_dir=temp_dir,
        spark_conf=spark_conf_more_ram)

hl.default_reference('GRCh38')

# DEFINE ALL PATHS
CHR_VAR = args.contig # 'chr20' # VARIABLE!!!

# INPUTS
AOU_VDS_PATH = args.input_aou_vds_url # 'gs://prod-drc-broad/v8/wgs/vds/aou_srwgs_short_variants_v8r1.vds/merged'

# OUTPUTS
vcf_url = args.output_aou_vcf_url
vcf_header_url = args.output_aou_vcf_header_url
rep_url = args.output_report_url

## IMPORT AOU DATA
VDS = hl.vds.read_vds(AOU_VDS_PATH)

## FILTER VDS TO CHROMOSOME (REGIONAL FILTERING)
VDS_FIL = hl.vds.filter_chromosomes(VDS, keep=CHR_VAR)

## DENSIFY THE VDS FILE TO MT
mt = hl.vds.to_dense_mt(VDS_FIL)

## FILTER BY THE NUMBER OF ALT ALLELES
ALT_MAX = 31 # VARIABLE!!!
mt = mt.filter_rows(hl.len(mt.alleles) < (ALT_MAX+2))

## ANNOTATE ENTRIES
mt = mt.annotate_entries(
    GT = hl.vds.lgt_to_gt(mt.LGT, mt.LA),
    AD = hl.vds.local_to_global(mt.LAD,
                                mt.LA,
                                n_alleles = hl.len(mt.alleles),
                                fill_value = 0,
                                number = 'R'))

mt = mt.annotate_entries(
    sumAD = hl.sum(mt.AD))

## ANNOTATE ROWS
mt = hl.variant_qc(mt)

mt = mt.annotate_rows(
    infor = hl.struct(AC = mt.variant_qc.AC[1:],
                      AF = mt.variant_qc.AF[1:],
                      AN = mt.variant_qc.AN,
                      homozygote_count = mt.variant_qc.homozygote_count))

mt = mt.annotate_rows(
    average_variant_sum_AD = hl.agg.sum(mt.sumAD)/hl.agg.count_where(mt.sumAD >= 0),
    maximum_variant_AC = hl.max(mt.infor.AC),
    defined_AD = hl.agg.count_where(mt.sumAD >= 0))

## DROP FIELDS
fields_to_drop = ['LAD', 'truth_sensitivity_snp_threshold', 'truth_sensitivity_indel_threshold', 'as_vets',
                  'LGT', 'LA', 'FT', 'PS', 'RGQ', 'LAD', 'sumAD', 'GQ', 'AD']

mt = mt.drop(*fields_to_drop)

## APPLY ROW FILTERS
# Filter the matrixtable
mt_fil = mt.filter_rows(
    ((mt.defined_AD >= 1) & (mt.average_variant_sum_AD < 12.0)) | # VARIABLE!!!
    (mt.maximum_variant_AC < 2) | # VARIABLE!!!
    ((mt.filters.contains('LowQual')) | # VARIABLE!!! Probably a checkbox?
     (mt.filters.contains('NO_HQ_GENOTYPES')) | # VARIABLE!!! Probably a checkbox?
     (mt.filters.contains('ExcessHet'))) | # VARIABLE!!! Probably a checkbox?
    (mt.variant_qc.call_rate < 0.9) |  # VARIABLE!!!
    (mt.variant_qc.gq_stats.mean < 30.0),  # VARIABLE!!!
    keep=False)

## ANNOTATE INFO FIELDS, REMOVE OTHER ETRANEOUS FIELDS
mt_fil = mt_fil.annotate_rows(
    info = hl.struct(GQ = mt_fil.variant_qc.gq_stats.mean,
                     AC = mt_fil.infor.AC,
                     AF = mt_fil.infor.AF,
                     AN = mt_fil.infor.AN,
                     HC = mt_fil.infor.homozygote_count,
                     AVSAD = mt_fil.average_variant_sum_AD))

fields_to_drop = ['filters', 'variant_qc', 'infor', 'maximum_variant_AC', 'defined_AD', 'average_variant_sum_AD']

mt_fil = mt_fil.drop(*fields_to_drop)

## CREATE AND STORE VCF HEADER
vcf_metadata = """##fileformat=VCFv4.2
##reference=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=HC,Number=R,Type=Integer,Description="Number of homozygotes per allele. One element per allele, including the reference.">
##INFO=<ID=AVSAD,Number=1,Type=Float,Description="Mean sum of allelic depts. Proxies DP.">
##INFO=<ID=GQ,Number=1,Type=Float,Description="Mean Genotype Quality">
"""
vcf_metadata += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t0000000000"

local_fname = "tmp_vcf_header.txt"
fp = open(local_fname, 'w')
fp.write(vcf_metadata)
fp.flush()
fp.close()

run_subprocess(f"gcloud storage cp {local_fname} {vcf_header_url}", f"Error copying header - {local_fname} to {vcf_header_url}")

## WRITE VCF
hl.export_vcf(mt_fil,
              vcf_url,
              metadata = hl.get_vcf_metadata(vcf_header_url),
              tabix = False)

## WRITE COST REPORT
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
cost_est = rig_cost * time_est

export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{vcf_url}\n" # @Franjo - in other scripts you use NAME_BASE instead of the path - does it matter?
    f"Time:\t{time_est}\n"
    f"Rigs:\t{rig_cost}\n"
    f"Cost:\t{cost_est}")

local_fname2 = "tmp_vcf_report.txt"
fp = open(local_fname2, 'w')
fp.write(export_report)
fp.flush()
fp.close()

run_subprocess(f"gcloud storage cp {local_fname2} {rep_url}", f"Error copying report - {local_fname2} to {rep_url}")