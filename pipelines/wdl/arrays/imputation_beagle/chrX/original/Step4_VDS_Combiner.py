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
parser.add_argument('--bed_url',                        help='GCS url to BED file for interval filtering')
# INPUTS
parser.add_argument('--input_aou_filtered_vds_url',     help="GCS url to the filtered AoU VDS input")
parser.add_argument('--input_anvil_filtered_vds_url',   help="GCS url to the filtered AnVIL CCDG VDS input")
# OUTPUTS
parser.add_argument('--output_combined_vds_url',        help="GCS url where the combined AoU-AnVIL VDS will be written")
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

# BOOT UP HAIL
temp_dir = f'{temp_bucket}/Step_4/temp'

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
BED_URL = args.bed_url

VDS_ANVIL_FIL_PATH = args.input_anvil_filtered_vds_url # f'{bucket}/Combined_Data/ANVIL_{NAME_BASE}.vds'
VDS_AOU_FIL_PATH = args.input_aou_filtered_vds_url # f'{bucket}/Combined_Data/AOU_{NAME_BASE}.vds'

COMBINED_PATH = args.output_combined_vds_url # f'{bucket}/Combined_Data/COMBINED_{NAME_BASE}.vds'

rep_url = args.output_report_url # f'{bucket}/Combined_Data/{CHR_VAR}.step4.report'

# Define the combiner config (adjust parameters as needed)
combiner_AOU = hl.vds.new_combiner(
    output_path=COMBINED_PATH,
    temp_path=temp_dir,
    vds_paths=[VDS_AOU_FIL_PATH, VDS_ANVIL_FIL_PATH],
    use_genome_default_intervals=True)

# Run the combiner
combiner_AOU.run()

# Cost
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
cost_est = rig_cost * time_est

# Write report
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{COMBINED_PATH}\n"
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
