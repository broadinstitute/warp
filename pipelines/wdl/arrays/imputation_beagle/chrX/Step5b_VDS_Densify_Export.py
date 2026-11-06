# SET UP THE ENVIRONMENT
from datetime import datetime
import os
import pandas as pd
import numpy as np
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--contig",                             help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',                        help='GCS bucket used for spark.local.dir, likely the cluster bucket')
parser.add_argument('--bed_url',                            help='GCS url to BED file for interval filtering')
# INPUTS
parser.add_argument('--input_filtered_vds_url',             help='GCS url to the combined AoU-AnVIL VDS input')
# OUTPUTS
parser.add_argument('--output_final_vcf_url',               help="GCS url where the combined AoU-AnVIL VCF will be written")
parser.add_argument('--output_final_vcf_header_url',        help="GCS url where the combined AoU-AnVIL VCF header will be written")
parser.add_argument('--output_report_url',                  help='GCS url where the report output will be copied')
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

# START THE TIME COUNTER
start = datetime.now()

# DEFINE BUCKET FOR THE WORKSPACE
temp_bucket = args.temp_bucket

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
rig_cost=15.06

# SET UP HAIL ENVIRONMENT (+ PLOTTING)
import hail as hl


# BOOT UP HAIL
temp_dir = f'{temp_bucket}/Step_5b/temp'

spark_conf_more_ram = dict()
spark_conf_more_ram["spark.executor.memory"] = "45g"
spark_conf_more_ram["spark.executor.cores"] = "8"
spark_conf_more_ram["spark.driver.memory"] = "300g"
spark_conf_more_ram["spark.driver.cores"] = "64"
spark_conf_more_ram["spark.local.dir"] = temp_dir

hl.init(idempotent=True,
        tmp_dir=temp_dir,
        spark_conf=spark_conf_more_ram)

hl.default_reference('GRCh38')

# DEFINE ALL PATHS
CHR_VAR = args.contig # 'chr20' # VARIABLE!!!
BED_URL = args.bed_url

FIL_VDS_PATH = args.input_filtered_vds_url

final_vcf_url = args.output_final_vcf_url
final_vcf_header_url = args.output_final_vcf_header_url
rep_url = args.output_report_url

# READ VDS files
VDS = hl.vds.read_vds(FIL_VDS_PATH)

# FILTER VDS TO CHROMOSOME (REGIONAL FILTERING OPTION 1)
VDS_X = hl.vds.filter_chromosomes(VDS, keep=CHR_VAR)
bed_intervals = hl.import_bed(BED_URL, reference_genome='GRCh38')
VDS_X = hl.vds.filter_intervals(VDS_X, bed_intervals, split_reference_blocks=True)


# DENSIFY THE VDS FILE TO MT
mt = hl.vds.to_dense_mt(VDS_X)

# FILTER BY THE NUMBER OF ALT ALLELES
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
fields_to_drop = ['LAD', 'LGT', 'LA', 'RGQ', 'LAD', 'sumAD', 'GQ', 'AD']

mt_fil = mt.drop(*fields_to_drop)

## ANNOTATE INFO FIELDS, REMOVE OTHER ETRANEOUS FIELDS
mt_fil = mt_fil.annotate_rows(
        info = hl.struct(GQ = mt_fil.variant_qc.gq_stats.mean,
                AC = mt_fil.infor.AC,
                AF = mt_fil.infor.AF,
                AN = mt_fil.infor.AN,
                HC = mt_fil.infor.homozygote_count,
                AVSAD = mt_fil.average_variant_sum_AD))

fields_to_drop = ['variant_qc', 'infor', 'maximum_variant_AC', 'defined_AD', 'average_variant_sum_AD']

mt_fil = mt_fil.drop(*fields_to_drop)

# Create and store VCF header
vcf_metadata = """##fileformat=VCFv4.2
##reference=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=homozygote_count,Number=R,Type=Integer,Description="Number of homozygotes per allele. One element per allele, including the reference.">
##INFO=<ID=average_variant_sum_AD,Number=1,Type=Float,Description="Mean sum of allelic depts. Proxies DP.">
##INFO=<ID=GQ,Number=1,Type=Float,Description="Mean Genotype Quality">
"""
# Fake a header line that ultimately gets ignored anyway.
vcf_metadata += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t0000000000"

# Save it to a GS url
local_fname = "tmp_vcf_header.txt"
fp = open(local_fname, 'w')
fp.write(vcf_metadata)
fp.flush()
fp.close()

# Copy the file to gs space
run_subprocess(f"gcloud storage cp {local_fname} {final_vcf_header_url}", f"Error copying header - {local_fname} to {final_vcf_header_url}")

# Write VCF
hl.export_vcf(mt_fil,
              final_vcf_url,
              metadata = hl.get_vcf_metadata(final_vcf_header_url),
              tabix = False)

# Cost
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
cost_est = rig_cost * time_est

# Write report
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{final_vcf_url}\n"
    f"Time:\t{time_est}\n"
    f"Rigs:\t{rig_cost}\n"
    f"Cost:\t{cost_est}")

# Save it to a GS url
local_fname2 = "tmp_vcf_report.txt"
fp = open(local_fname2, 'w')
fp.write(export_report)
fp.flush()
fp.close()

# Copy the file to gs space
run_subprocess(f"gcloud storage cp {local_fname2} {rep_url}", f"Error copying report - {local_fname2} to {rep_url}")
