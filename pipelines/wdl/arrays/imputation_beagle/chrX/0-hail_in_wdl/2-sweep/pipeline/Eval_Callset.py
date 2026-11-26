"""Evaluation callset: unfiltered, unphased genotypes of a few AoU samples.

For the samples that have long-read truth, write a small VCF with their genotypes
at EVERY AoU site in the region where at least one of them carries a non-reference
allele (sites with more than ALT_MAX alt alleles excepted). No site filter is
applied, so one vcfdist run against the truth scores every candidate site, and
filter_sweep.py can then score any filter set offline.

Genotypes are made the way Step 4 makes them: global GT, haploid calls recoded as
homozygous diploid, densified reference blocks as 0/0. They are unphased; run
VcfdistEvaluation.wdl on this VCF with do_naively_phase = true.
"""
import argparse
import subprocess
import sys
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--contig',                     required=True, help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',                required=True, help='GCS path used for Hail temp files, likely the cluster temp bucket')
parser.add_argument('--bed_url',                    required=True, help='GCS url to BED file for interval filtering')
parser.add_argument('--samples',                    required=True, help='comma-separated AoU sample IDs')
parser.add_argument('--ref_pad_bp',                 type=int, default=1000000)
parser.add_argument('--executor_memory',            default='22g')
parser.add_argument('--executor_cores',             default='4')
parser.add_argument('--rig_cost',                   type=float, default=0.0)
# INPUTS
parser.add_argument('--input_aou_vds_url',          required=True, help='GCS url to the AoU VDS')
# OUTPUTS
parser.add_argument('--output_vcf_url',             required=True, help='GCS url of the evaluation VCF (.vcf.bgz)')
parser.add_argument('--output_report_url',          required=True, help='GCS url where the report output will be copied')
args = parser.parse_args()

ALT_MAX = 31
SAMPLES = sorted({s for s in args.samples.split(',') if s})
if not SAMPLES:
    sys.exit("--samples is empty")


def run_subprocess(cmd, errorMessage):
    if isinstance(cmd, list):
        cmd = ' '.join(cmd)
    try:
        return subprocess.check_output(cmd, shell=True, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        print(errorMessage)
        print("Exited with " + str(e.returncode) + "-" + e.output)
        exit(1)


start = datetime.now()

import hail as hl  # noqa: E402

temp_dir = f'{args.temp_bucket}/Eval_Callset/temp'
hl.init(idempotent=True,
        tmp_dir=temp_dir,
        spark_conf={"spark.executor.memory": args.executor_memory,
                    "spark.executor.cores": args.executor_cores,
                    "spark.driver.memory": "60g",
                    "spark.driver.cores": "32",
                    "spark.local.dir": temp_dir})
hl.default_reference('GRCh38')

CHR_VAR = args.contig


# ---------------------------------------------------------------- shared helpers
def load_intervals(bed_url, contig):
    """BED intervals on one contig, as Python hl.Interval objects sorted by start."""
    bed = hl.import_bed(bed_url, reference_genome='GRCh38')
    ivs = bed.filter(bed.interval.start.contig == contig).interval.collect()
    if not ivs:
        sys.exit(f"No intervals on {contig} in {bed_url}")
    return sorted(ivs, key=lambda iv: iv.start.position)


def read_region(vds_path, contig, intervals, ref_pad_bp, label):
    """Read a VDS restricted to the intervals, pruning partitions on both halves.

    Variant data: rows inside the intervals. Reference data: rows that start inside
    the intervals or up to `pad` bp before them, so that blocks which begin earlier
    and reach into an interval are still there for densifying. `pad` is the stored
    maximum reference-block length when the VDS has one, otherwise ref_pad_bp.
    hl.vds.filter_intervals is deliberately not used: without a stored maximum it
    leaves the whole contig's reference data in the plan (a full pass).
    """
    vds = hl.vds.read_vds(vds_path)
    rd, vd = vds.reference_data, vds.variant_data
    fd = hl.vds.VariantDataset.ref_block_max_length_field
    if fd in rd.globals:
        pad = int(hl.eval(rd.index_globals()[fd])) - 1
        print(f"{label}: stored maximum reference-block length -> padding {pad} bp")
    else:
        pad = ref_pad_bp
        print(f"{label}: no stored maximum reference-block length -> padding {pad} bp (--ref_pad_bp)")
    padded = [hl.Interval(hl.Locus(contig, max(1, iv.start.position - pad), reference_genome='GRCh38'),
                          iv.end, includes_start=True, includes_end=iv.includes_end)
              for iv in intervals]
    rd = hl.filter_intervals(rd, padded)
    vd = hl.filter_intervals(vd, intervals)
    try:
        print(f"{label}: {vd.n_partitions()} variant partitions, {rd.n_partitions()} reference partitions after interval filter")
    except Exception as exc:  # informational only
        print(f"{label}: could not count partitions ({exc})")
    return rd, vd


intervals = load_intervals(args.bed_url, CHR_VAR)
rd, vd = read_region(args.input_aou_vds_url, CHR_VAR, intervals, args.ref_pad_bp, 'AoU')

wanted = hl.literal(set(SAMPLES))
vd = vd.filter_cols(wanted.contains(vd.s))
rd = rd.filter_cols(wanted.contains(rd.s))
found = vd.s.collect()
missing = sorted(set(SAMPLES) - set(found))
if missing:
    sys.exit(f"{len(missing)} of {len(SAMPLES)} samples are not in the VDS: {missing[:10]}")
print(f"{len(found)} evaluation samples found")

vd = vd.filter_rows(hl.len(vd.alleles) < (ALT_MAX + 2))
# Global GT with haploid calls recoded as homozygous diploid, as in Step 4
gt = hl.vds.lgt_to_gt(vd.LGT, vd.LA)
vd = vd.select_entries(GT=hl.if_else(gt.ploidy == 1, hl.call(gt[0], gt[0]), gt))
# Keep the sites where an evaluation sample carries a non-reference allele
vd = vd.filter_rows(hl.agg.any(vd.GT.is_non_ref()))
vd = vd.select_rows()
rd = rd.select_entries(*[f for f in ('END', 'LEN') if f in rd.entry])

mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(rd, vd))
mt = mt.select_entries('GT')
hl.export_vcf(mt, args.output_vcf_url)

# ----------------------------------------------------------------------- report
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{args.output_vcf_url}\n"
    f"Time:\t{time_est}\n"
    f"Rigs:\t{args.rig_cost}\n"
    f"Cost:\t{args.rig_cost * time_est}")
local_fname = "tmp_step_report.txt"
with open(local_fname, 'w') as fp:
    fp.write(export_report)
run_subprocess(f"gcloud storage cp {local_fname} {args.output_report_url}",
               f"Error copying report - {local_fname} to {args.output_report_url}")
