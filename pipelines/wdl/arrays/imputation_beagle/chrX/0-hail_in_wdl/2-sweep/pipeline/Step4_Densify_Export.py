"""Step 4: densify a VDS and export it as the final VCF (was Step 5b).

The input is the combined AoU-AnVIL VDS from Step 3. It can equally be Step 3's
filtered AoU-only VDS, which gives an AoU-only VCF at exactly the same kept sites
(the comparison arm for "what does adding AnVIL buy").

The VCF has the same records, INFO fields (GQ, AC, AF, AN, HC, AVSAD) and GT
values as the original Step5b_VDS_Densify_Export: haploid calls are recoded as
homozygous diploid before any counting, and sites with more than ALT_MAX alt
alleles are dropped.

One deliberate difference from the original records, on by default
(--mask_xy_nonpar_hets true, chrX only): a heterozygous call in an XY sample in
chrX non-PAR is written as missing (./.), so Beagle imputes it. A male has one
X there, so such a call cannot be right; AnVIL males are diploid-encoded and
typically have 5 to 10% of their non-PAR non-reference calls heterozygous. AoU
males are haploid-encoded and are unaffected. AC / AF / AN / HC are counted after
masking, so they describe the genotypes in the file. With the option false the
records are those of the original.

What is otherwise different from the original, without changing the records:
  * No interval filter: the combined VDS already contains only the kept sites, and
    filtering it again made Hail walk the reference data to split blocks.
  * No naive_coalesce(16): that ran the whole densify-and-export on 16 tasks
    regardless of cluster size. The export now uses the combined VDS's own
    partitions (the tabix index is still written during the export).
  * AVSAD comes from the sparse variant data (only variant calls carry AD), and
    the allele-depth array is no longer rebuilt per densified entry.
  * GT conversion and haploid recoding happen on the sparse variant data, so the
    dense matrix carries two entry fields (GT, GQ).
  * hl.variant_qc is replaced by the two aggregations whose results are written
    (call_stats and GQ stats); these are the same aggregators variant_qc uses.
"""
import argparse
import re
import subprocess
import sys
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--contig',                             required=True, help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',                        required=True, help='GCS path used for Hail temp files, likely the cluster temp bucket')
parser.add_argument('--executor_memory',                    default='45g')
parser.add_argument('--executor_cores',                     default='8')
parser.add_argument('--rig_cost',                           type=float, default=0.0, help='cluster cost per hour, for the report only')
parser.add_argument('--mask_xy_nonpar_hets',                default='true', choices=['true', 'false'], help='chrX only: write heterozygous calls of XY samples in non-PAR as missing')
parser.add_argument('--aou_sex_tsv_url',                    default='', help='two-column (sample, sex) table for AoU; needed for masking')
parser.add_argument('--anvil_sex_tsv_url',                  default='', help='two-column (sample, sex) table for AnVIL, with the sample IDs as in the AnVIL VDS; needed for masking')
parser.add_argument('--anvil_sample_suffix',                default='_C', help='suffix Step 3 appends to AnVIL sample IDs')
# INPUTS
parser.add_argument('--input_combined_vds_url',             required=True, help='combined AoU-AnVIL VDS (or the filtered AoU-only VDS), from Step 3')
# OUTPUTS
parser.add_argument('--output_final_vcf_url',               required=True, help='GCS url where the combined AoU-AnVIL VCF will be written')
parser.add_argument('--output_final_vcf_header_url',        required=True, help='GCS url where the combined AoU-AnVIL VCF header will be written')
parser.add_argument('--output_report_url',                  required=True, help='GCS url where the report output will be copied')
args = parser.parse_args()

ALT_MAX = 31
NONPAR_START, NONPAR_END = 2781480, 155701382      # GRCh38 chrX
XY_LABELS = {'XY', 'M', 'MALE'}
MASK_XY_HETS = args.mask_xy_nonpar_hets == 'true' and args.contig == 'chrX'
if MASK_XY_HETS and not (args.aou_sex_tsv_url or args.anvil_sex_tsv_url):
    sys.exit("--mask_xy_nonpar_hets true needs --aou_sex_tsv_url and/or --anvil_sex_tsv_url")


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

temp_dir = f'{args.temp_bucket}/Step_4/temp'
hl.init(idempotent=True,
        tmp_dir=temp_dir,
        spark_conf={"spark.executor.memory": args.executor_memory,
                    "spark.executor.cores": args.executor_cores,
                    "spark.driver.memory": "300g",
                    "spark.driver.cores": "64",
                    "spark.local.dir": temp_dir})
hl.default_reference('GRCh38')

CHR_VAR = args.contig

vds = hl.vds.read_vds(args.input_combined_vds_url)
rd, vd = vds.reference_data, vds.variant_data
try:
    print(f"combined VDS: {vd.n_partitions()} variant partitions, {rd.n_partitions()} reference partitions")
except Exception as exc:  # informational only
    print(f"could not count partitions ({exc})")

# ------------------------------------------------------------- XY samples
def load_xy_ids(sex_tsv_url, suffix):
    ids = set()
    if sex_tsv_url:
        with hl.hadoop_open(sex_tsv_url, 'r') as fh:
            for line in fh:
                parts = re.split(r'[\t, ]+', line.strip())
                if len(parts) >= 2 and parts[1].upper() in XY_LABELS:
                    ids.add(parts[0] + suffix)
    return ids


if MASK_XY_HETS:
    aou_xy = load_xy_ids(args.aou_sex_tsv_url, '')
    anvil_xy = load_xy_ids(args.anvil_sex_tsv_url, args.anvil_sample_suffix)
    xy_ht = hl.Table.parallelize([{'s': s} for s in sorted(aou_xy | anvil_xy)],
                                 schema=hl.tstruct(s=hl.tstr), key='s')
    vd = vd.annotate_cols(_is_xy=hl.is_defined(xy_ht[vd.s]))
    n_xy = vd.aggregate_cols(hl.agg.count_where(vd._is_xy))
    print(f"XY labels: {len(aou_xy)} AoU, {len(anvil_xy)} AnVIL (with suffix '{args.anvil_sample_suffix}'); "
          f"{n_xy} of {vd.count_cols()} samples in this VDS matched")
    if n_xy == 0:
        sys.exit("No sample in the VDS matched an XY label; check the sex tables and --anvil_sample_suffix")

# ------------------------------------------------------ sparse (variant data)
vd = vd.filter_rows(hl.len(vd.alleles) < (ALT_MAX + 2))

# Sum of the global AD array == sum of the local one (local_to_global only adds zeros)
sum_ad = hl.or_missing(hl.is_defined(vd.LA), hl.sum(vd.LAD))
vd = vd.annotate_rows(_sum_ad=hl.agg.sum(sum_ad), _n_ad=hl.agg.count_where(sum_ad >= 0))

# Global GT, with haploid calls recoded as homozygous diploid (1 -> 1/1). Densified
# reference entries are filled with diploid 0/0, so every sample ends up with one
# ploidy, as Beagle requires.
gt = hl.vds.lgt_to_gt(vd.LGT, vd.LA)
gt = hl.if_else(gt.ploidy == 1, hl.call(gt[0], gt[0]), gt)
if MASK_XY_HETS:
    # A het call in an XY sample in non-PAR cannot be right: write it as missing. The entry
    # stays a variant-data entry, so densifying leaves it missing instead of filling 0/0.
    in_nonpar = (vd.locus.position >= NONPAR_START) & (vd.locus.position <= NONPAR_END)
    gt = hl.or_missing(~(in_nonpar & vd._is_xy & gt.is_het()), gt)
vd = vd.select_entries(GT=gt, GQ=vd.GQ)
rd = rd.select_entries(*([f for f in ('END', 'LEN') if f in rd.entry] + ['GQ']))

# ------------------------------------------------------------------- densify
mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(rd, vd))
mt = mt.annotate_rows(_cs=hl.agg.call_stats(mt.GT, mt.alleles), _gq=hl.agg.stats(mt.GQ))
mt = mt.select_rows(
    info=hl.struct(GQ=mt._gq.mean,
                   AC=mt._cs.AC[1:],
                   AF=mt._cs.AF[1:],
                   AN=mt._cs.AN,
                   HC=mt._cs.homozygote_count,
                   AVSAD=mt._sum_ad / mt._n_ad))
mt = mt.select_entries('GT')
mt = mt.select_cols()

# ---------------------------------------------------------------- VCF header
# Unchanged from the original Step 5b.
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

local_fname = "tmp_vcf_header.txt"
with open(local_fname, 'w') as fp:
    fp.write(vcf_metadata)
run_subprocess(f"gcloud storage cp {local_fname} {args.output_final_vcf_header_url}",
               f"Error copying header - {local_fname} to {args.output_final_vcf_header_url}")

# Write VCF, building the tabix index during the export (<final_vcf_url>.tbi)
hl.export_vcf(mt,
              args.output_final_vcf_url,
              metadata=hl.get_vcf_metadata(args.output_final_vcf_header_url),
              tabix=True)

# ----------------------------------------------------------------------- report
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{args.output_final_vcf_url}\n"
    f"Time:\t{time_est}\n"
    f"Rigs:\t{args.rig_cost}\n"
    f"Cost:\t{args.rig_cost * time_est}")
local_fname2 = "tmp_step_report.txt"
with open(local_fname2, 'w') as fp:
    fp.write(export_report)
run_subprocess(f"gcloud storage cp {local_fname2} {args.output_report_url}",
               f"Error copying report - {local_fname2} to {args.output_report_url}")
