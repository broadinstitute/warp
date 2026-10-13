# Given a split dense Hail MT (post sc4_export_dense_bed_MTs_to_split_dense_MTs.py), export it as a VCF
# per callset/chromosome. This is the first time this pipeline produces a VCF output -- there is no
# equivalent VDS-input script to mirror.
#
# PROTOTYPE STATUS: the header/metadata builder and the core per-MT export function are implemented and
# directly testable against an existing MT (e.g. the sc1_create_basis_mt_from_vcf.py chr21 test output),
# without needing SC4 to exist yet. Still open before this is a real, wired-up export script:
#   - Whether the LowQual/ExcessHet/NO_HQ_GENOTYPES FILTER definitions below are ever actually needed
#     (SC4 currently drops any row with a non-empty `filters` field before export)
#   - FORMAT field coverage beyond GT/AD -- need to check the real phased VCF's FORMAT columns
#     (bcftools view -h) for anything like a phase-set field that would need its own header entry
#   - Whether `##reference=` needs to be supplied explicitly or Hail fills it in automatically
#   - Reading SC4's actual split-dense-MT json manifest and looping per bed_name (not built yet --
#     SC4's output shape isn't settled, so this would be guessing at an interface)
#   - Per-chromosome output naming, indexing (tabix/bcftools, not Hail's own tabix=True), CLI/WDL wiring
#
# See all_of_us/phasing/AoU_VCF_v9.py and all_of_us/phasing/remove_phased_samples_X.py
# (build_default_metadata) for the reference conventions this is based on.


import hail as hl


def require_gs_path(path: str, label: str) -> None:
    """
    Raise if the given path is not a gs:// bucket path. Same rationale as
    sc1_create_basis_mt_from_vcf.py: a VCF export has to be written directly to a bucket.

    Parameters:
    path (str): The path to check.
    label (str): A human-readable label for the path, used in the error message.
    """
    if not path.startswith("gs://"):
        raise ValueError(f"{label} must be a gs:// bucket path -- got: {path}")


def check_output_does_not_exist(path: str, overwrite_existing_outputs: bool) -> None:
    """
    Guard against silently clobbering an existing output. Same convention as
    sc1_create_basis_mt_from_vcf.py.

    Parameters:
    path (str): The gs:// path to check.
    overwrite_existing_outputs (bool): If True, skip the check.

    Raises:
    FileExistsError: If `path` already exists and overwrite_existing_outputs is False.
    """
    if overwrite_existing_outputs:
        return
    if hl.hadoop_exists(path):
        raise FileExistsError(
            f"Output already exists at {path} and overwrite_existing_outputs is False. "
            f"Pick a new output path, or pass overwrite_existing_outputs=True if you really mean to "
            f"replace it."
        )


def build_phased_vcf_metadata() -> dict:
    """
    VCF header metadata for the phased small-callset export.

    AC/AF/AN/HC are recomputed by this pipeline (post BED-filter, post split_multi_hts). GQ and AVSAD
    are carried forward unchanged from the original (pre-phasing) INFO field -- matching the convention
    established in all_of_us/phasing/AoU_VCF_v9.py and all_of_us/phasing/remove_phased_samples_X.py.

    FORMAT only declares GT -- confirmed via a real exported header that this phased VCF carries no other
    FORMAT fields (no AD, no phase-set field) worth declaring. Revisit if that ever changes.

    Returns:
    dict: suitable for hl.export_vcf(..., metadata=build_phased_vcf_metadata()).
    """
    return {
        "info": {
            "AC": {
                "Number": "A",
                "Type": "Integer",
                "Description": "Allele count in genotypes, for each ALT allele, in the same order as listed",
            },
            "AF": {
                "Number": "A",
                "Type": "Float",
                "Description": "Allele Frequency, for each ALT allele, in the same order as listed",
            },
            "AN": {
                "Number": "1",
                "Type": "Integer",
                "Description": "Total number of alleles in called genotypes",
            },
            "HC": {
                "Number": "R",
                "Type": "Integer",
                "Description": "Number of homozygotes per allele. One element per allele, including the reference.",
            },
            "AVSAD": {
                "Number": "1",
                "Type": "Float",
                "Description": "Mean sum of allelic depths. Proxies DP. Carried forward from the "
                               "pre-phasing callset, not recomputed by this pipeline.",
            },
            "GQ": {
                "Number": "1",
                "Type": "Float",
                "Description": "Mean Genotype Quality. Carried forward from the pre-phasing callset, "
                               "not recomputed by this pipeline.",
            },
        },
        "filter": {
            "LowQual": {"Description": "Low quality score"},
            "ExcessHet": {"Description": "Excess heterozygotes"},
            "NO_HQ_GENOTYPES": {"Description": "No high-quality genotypes"},
        },
        "format": {
            "GT": {
                "Number": "1",
                "Type": "String",
                "Description": "Genotype",
            },
        },
    }


def export_mt_to_vcf(mt: hl.MatrixTable, vcf_url: str, overwrite_existing_outputs: bool = False) -> None:
    """
    Export a Hail MatrixTable to a VCF, using this pipeline's header/metadata convention. Indexing is
    intentionally NOT done here (Hail's own tabix=True) -- per guidance, index separately with a cheaper
    non-Hail tool (tabix/bcftools) as its own task, not as part of this Dataproc/Hail step.

    Parameters:
    mt (hl.MatrixTable): The MatrixTable to export. Expected to already be BED-filtered, phase-preserved,
                        and have info.AC/AF/AN/HC recomputed (e.g. the output of
                        sc1_create_basis_mt_from_vcf.py, and eventually sc4's split/recompute step).
    vcf_url (str): gs:// path to write the VCF to (e.g. ending in .vcf.bgz).
    overwrite_existing_outputs (bool): If False (default), refuse to run if vcf_url already exists.

    Returns:
    None
    """
    require_gs_path(vcf_url, "vcf_url")
    check_output_does_not_exist(vcf_url, overwrite_existing_outputs)
    hl.export_vcf(mt, vcf_url, metadata=build_phased_vcf_metadata(), tabix=False)
