# SET UP THE ENVIRONMENT
from datetime import datetime
import os
import pandas as pd
import numpy as np
import argparse
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--contig',                 help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',            help='GCS bucket used for spark.local.dir, likely the cluster bucket')
# INPUTS
parser.add_argument('--input_anvil_vds_url',    help='GCS url to the AnVIL CCDG-only VDS file to process (e.g. gs://bucket/path/to/file.vds)')
# OUTPUTS
parser.add_argument('--output_anvil_mt_url',    help='GCS url where the AnVIL CCDG-only MT file will be copied')
parser.add_argument('--output_report_url',      help='GCS url where the report output will be copied')
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

# DEFINE BUCKETS
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
rig_cost=7.80

# SET UP HAIL ENVIRONMENT (+ PLOTTING)
import hail as hl

# BOOT UP HAIL
temp_dir = f'{temp_bucket}/Step_2/temp'

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
ANVIL_VDS_PATH = args.input_anvil_vds_url # 'gs://fc-secure-759ddb91-68ac-4234-8fab-9c9325c3484f/CCDG_consent.vds'

# OUTPUTS
rep_url = args.output_report_url
MT_ANVIL_PATH = args.output_anvil_mt_url

# IMPORT ANVIL DATA
VDS = hl.vds.read_vds(ANVIL_VDS_PATH)

# FILTER VDS TO CHROMOSOME (REGIONAL FILTERING OPTION 1)
VDS_FIL = hl.vds.filter_chromosomes(VDS, keep=CHR_VAR)

# DENSIFY THE VDS FILE TO MT
mt = hl.vds.to_dense_mt(VDS_FIL)

# FILTER BY THE NUMBER OF ALT ALLELES
ALT_MAX = 31 # VARIABLE!!!
mt = mt.filter_rows(hl.len(mt.alleles) < (ALT_MAX+2))

# ANNOTATE ENTRIES
mt = mt.annotate_entries(
        GT = hl.vds.lgt_to_gt(mt.LGT, mt.LA),
        AD = hl.vds.local_to_global(mt.LAD,
                                mt.LA,
                                n_alleles = hl.len(mt.alleles),
                                fill_value = 0,
                                number = 'R'))

mt = mt.annotate_entries(
    sumAD = hl.sum(mt.AD))

NO_HQ_GENOTYPES_COND = ((mt.GQ>=20) & (mt.DP>=10) & 
                           ((mt.GT.is_het_ref() & ((mt.AD[1] / hl.sum(mt.AD))>=0.2)) |
                           (mt.GT.is_het_non_ref() & ((mt.AD[0] / hl.sum(mt.AD))>=0.2) & ((mt.AD[1] / hl.sum(mt.AD))>=0.2)) |
                           (~mt.GT.is_het_ref() & ~mt.GT.is_het_non_ref())))

mt=mt.annotate_entries(Pass_NO_HQ_GENOTYPES=hl.if_else(NO_HQ_GENOTYPES_COND, 1, 0))

# ANNOTATE ROWS
mt = hl.variant_qc(mt)

mt = mt.annotate_rows(
    info = hl.struct(AC = mt.variant_qc.AC[1:],
                     AF = mt.variant_qc.AF[1:],
                     AN = mt.variant_qc.AN,
                     homozygote_count = mt.variant_qc.homozygote_count))

mt = mt.annotate_rows(
    NO_HQ_GENOTYPES=hl.if_else(hl.agg.sum(mt.Pass_NO_HQ_GENOTYPES)>0, 0, 1),
    EXCESS_HET=hl.if_else(mt.variant_qc.p_value_excess_het>3.4e-06, 0, 1),
    average_variant_sum_AD = hl.agg.sum(mt.sumAD)/hl.agg.count_where(mt.sumAD >= 0),
    maximum_variant_AC = hl.max(mt.info.AC),
    defined_AD = hl.agg.count_where(mt.sumAD >= 0))

# DROP FIELDS
fields_to_drop = ['LGT', 'LAD', 'LA', 'LPGT', 'LPL', 'MIN_DP', 'SB', 'PID', 'gvcf_info']

mt = mt.drop(*fields_to_drop)

# APPLY ROW FILTERS
# Filter the matrixtable
mt_fil = mt.filter_rows(
    ((mt.defined_AD >= 1) & (mt.average_variant_sum_AD < 12.0)) | # VARIABLE!!!
    #(mt.maximum_variant_AC < 2) | # VARIABLE!!!
    hl.coalesce(mt.variant_qc.p_value_excess_het < 3.4e-6, False) | # VARIABLE!!! Probably a checkbox?
    hl.coalesce(hl.agg.sum(mt.Pass_NO_HQ_GENOTYPES) == 0, False) | # VARIABLE!!! Probably a checkbox?
    (mt.variant_qc.call_rate < 0.9) |  # VARIABLE!!!
    (mt.variant_qc.gq_stats.mean < 30.0),  # VARIABLE!!!
    keep=True)

# Write mt
mt_fil.write(MT_ANVIL_PATH, overwrite = True)

# Cost
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
cost_est = rig_cost * time_est

# Write report
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{MT_ANVIL_PATH}\n"
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
