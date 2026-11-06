# SET UP THE ENVIRONMENT
from datetime import datetime
import os
import pandas as pd
import numpy as np
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--contig',                         help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',                    help='GCS bucket used for spark.local.dir, likely the cluster bucket')
# INPUTS
parser.add_argument('--input_anvil_filtered_vds_url',   help='GCS path to the filtered AnVIL CCDG-only VDS file to process (e.g. gs://bucket/path/to/file.vds)')
# OUTPUTS
parser.add_argument('--output_anvil_renamed_vds_url',   help="GCS url where the renamed AnVIL CCDG VDS will be written")
parser.add_argument('--output_report_url',              help='GCS url where the report output will be copied')
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
rig_cost=7.80

# SET UP HAIL ENVIRONMENT (+ PLOTTING)
import hail as hl

temp_dir = f'{temp_bucket}/Step_3b/temp'

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
CHR_VAR = args.contig

VDS_ANVIL_FILTERED_PATH = args.input_anvil_filtered_vds_url
VDS_ANVIL_RENAMED_PATH = args.output_anvil_renamed_vds_url

rep_url = args.output_report_url

# READ VDS files
VDS_ANVIL = hl.vds.read_vds(VDS_ANVIL_FILTERED_PATH)

# FILTER VDS TO CHROMOSOME (REGIONAL FILTERING OPTION 1)
VDS_ANVIL_X = hl.vds.filter_chromosomes(VDS_ANVIL, keep=CHR_VAR)

# Drop fields that are not present in both ANVIL data
VDS_FIL_ANVIL_VAR = VDS_ANVIL_X.variant_data
VDS_FIL_ANVIL_REF = VDS_ANVIL_X.reference_data

# Rename ANVIL Samples
VDS_FIL_ANVIL_VAR = VDS_FIL_ANVIL_VAR.key_cols_by(s = VDS_FIL_ANVIL_VAR.s.replace('$', '_C'))
VDS_FIL_ANVIL_REF = VDS_FIL_ANVIL_REF.key_cols_by(s = VDS_FIL_ANVIL_REF.s.replace('$', '_C'))

# Recombine variant and reference data
VDS_FILCOM_ANVIL = hl.vds.VariantDataset(VDS_FIL_ANVIL_REF, VDS_FIL_ANVIL_VAR)

# Write VDS file for ANVIL
VDS_FILCOM_ANVIL.write(VDS_ANVIL_RENAMED_PATH, overwrite=True)

# Cost
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
cost_est = rig_cost * time_est

# Write report
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{VDS_ANVIL_RENAMED_PATH}\n"
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

