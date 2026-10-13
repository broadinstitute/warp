# Given one or more phased, per-chromosome VCFs and a list of bed files, create a matrixTable of
# filtered variants to regions in these bed files.
# All bed files referenced in this python script are UCSC bed files (as opposed to PLINK bed files).
#
# This is the phased-VCF-input counterpart to sc1_create_basis_mt_from_vds.py. The VDS version reads a
# sparse Variant Dataset and has to reconstruct dense GT/AD from the local-allele (LGT/LA) encoding
# before it can do anything else. A phased VCF is already a dense, per-sample-called MatrixTable, so
# that reconstruction step goes away entirely; what's left is the format-agnostic part of the original
# script (BED annotate/filter, variant_qc, INFO recompute).
#
# PROTOTYPE STATUS: this is written for interactive testing on the hail fat box against a small region /
# single chromosome before it gets wired into a WDL. See sc1_create_basis_mt_from_vds.py for the
# production VDS-input sibling and the small_callsets README for pipeline context.


import argparse
import hail as hl
from typing import Dict, List, Optional, Union
import json
import subprocess
import re


def hail_init(executor_memory: str, executor_cores: str, driver_memory: str, reference_genome: str) -> None:
    """
    Initialize Hail with specific Spark configuration settings.

    Parameters:
    executor_memory (str): Memory assigned to each Spark executor, e.g., "15g".
    executor_cores (str): Number of cores assigned to each Spark executor, e.g., "8".
    driver_memory (str): Memory assigned to the Spark driver, e.g., "200".
    reference_genome (str): Reference genome identifier, e.g., "GRCh38".
    """
    spark_conf = {
        "spark.executor.memory": executor_memory,
        "spark.executor.cores": executor_cores,
        "spark.driver.memory": driver_memory
    }
    hl.init(default_reference=reference_genome, idempotent=True, spark_conf=spark_conf,
            quiet=True, skip_logging_configuration=True)


def create_bed_dict(exome_bed_path: str, clinvar_bed_path: str, acaf_threshold_bed_path: str) -> Dict[str, str]:
    """
    Create a dictionary mapping UCSC BED file types to their file paths.

    Parameters:
    exome_bed_path (str): Path to the exome UCSC BED file.
    clinvar_bed_path (str): Path to the ClinVar BED file.
    acaf_threshold_bed_path (str): Path to the aCAF threshold UCSC BED file.

    Returns:
    Dict[str, str]: A dictionary mapping UCSC BED file types to paths.
    """
    return {
        "exome": exome_bed_path,
        "clinvar": clinvar_bed_path,
        "acaf_threshold": acaf_threshold_bed_path
    }


def parse_vcf_paths(vcf_url: str) -> Union[List[str], str]:
    """
    Parse the --vcf_url argument into what hl.import_vcf expects. Supports either a single path/glob
    (e.g. a wildcard over per-chromosome files) or an explicit comma/whitespace-delimited list of paths
    (e.g. "gs://.../chr21.vcf.bgz,gs://.../chr22.vcf.bgz") for combining a few chromosomes in one run.

    Parameters:
    vcf_url (str): One path/glob, or a delimited list of paths.

    Returns:
    List[str] or str: A single path/glob string, or a list of paths, suitable for hl.import_vcf.
    """
    parts = [p for p in re.split(r'[,\s]+', vcf_url.strip()) if p]
    return parts[0] if len(parts) == 1 else parts


def import_phased_vcf(vcf_url: str, reference_genome: str, min_partitions: int) -> hl.MatrixTable:
    """
    Import one or more phased VCFs into a Hail MatrixTable.

    Parameters:
    vcf_url (str): Path/glob or delimited list of paths to the input phased VCF(s), bgzipped and tabix'd.
    reference_genome (str): Reference genome identifier, e.g., "GRCh38".
    min_partitions (int): Minimum number of partitions to import with, for parallelization.

    Returns:
    hl.MatrixTable: The imported MatrixTable. GT calls parsed from a "|"-delimited VCF genotype come in
                    already phased (hl.Call.phased == True) -- no extra work is needed at import time to
                    preserve phase.
    """
    paths = parse_vcf_paths(vcf_url)
    mt = hl.import_vcf(
        paths,
        reference_genome=reference_genome,
        force_bgz=True,
        array_elements_required=False,
        min_partitions=min_partitions,
    )
    return mt


def require_gs_path(path: str, label: str) -> None:
    """
    Raise if the given path is not a gs:// bucket path.

    A Hail MatrixTable is a directory of many partition files, not a single blob -- it must be written
    directly to a bucket by Hail. There's no way to write one to local (VM) disk and `gsutil cp` it to a
    bucket afterward, so the output location has to be a bucket path up front. Hail also needs to read
    VCF input directly from a bucket path rather than a local copy.

    Parameters:
    path (str): The path to check.
    label (str): A human-readable label for the path, used in the error message.
    """
    if not path.startswith("gs://"):
        raise ValueError(f"{label} must be a gs:// bucket path -- got: {path}")


def subset_to_samples(mt: hl.MatrixTable, sample_ids: Optional[str]) -> hl.MatrixTable:
    """
    Optionally subset the MatrixTable to a small, explicit list of sample/person IDs, e.g. for a fast
    fat-box test on 5 samples instead of the full callset. This mirrors the keep/remove-list idiom used
    in all_of_us/phasing/remove_phased_samples_X.py (import a list of IDs, then filter_cols on membership
    against mt.s) -- inverted here to keep instead of remove. For a longer list than a handful of IDs,
    switch to that file's hl.import_table(...)-based pattern instead of an inline literal.

    Parameters:
    mt (hl.MatrixTable): The input MatrixTable.
    sample_ids (Optional[str]): Comma/space-delimited list of sample IDs to keep, or None/empty to skip.

    Returns:
    hl.MatrixTable: The optionally sample-subsetted MatrixTable.
    """
    if sample_ids:
        sample_ids_list = [s for s in re.split(r'[,\s]+', sample_ids.strip()) if s]
        mt = mt.filter_cols(hl.literal(set(sample_ids_list)).contains(mt.s))
        print(f"Subsetting to {len(sample_ids_list)} sample(s): {sample_ids_list}")
    else:
        print("No sample subsetting applied.")
    return mt


def validate_input_chromosomes(test_chromosomes: list, reference_chromosomes: list) -> bool:
    """
    Validate input chromosomes against a reference list.

    Parameters:
    test_chromosomes (List): List of input chromosomes to be validated.
    reference_chromosomes (List): List of reference chromosomes to compare against.

    Returns:
    bool: Returns True if all chromosomes in 'test_chromosomes' are valid and present in 'reference_chromosomes'.
          Returns False if 'test_chromosomes' contains any invalid chromosomes not present in 'reference_chromosomes'
                        and exits analysis.
    """
    for value in test_chromosomes:
        if value not in reference_chromosomes:
            print("Input test chromosomes '{}' contains invalid chromosomes. Exiting analysis.".format(value))
            return False
    return True


def perform_chr_filtering_if_needed(mt: hl.MatrixTable, is_test: bool, test_chrs: str) -> hl.MatrixTable:
    """
    Optionally filter the MatrixTable down to specific chromosomes. Since phased input is already split
    per chromosome, this is mostly useful when --vcf_url was given a list/glob spanning more than one
    chromosome and you want to smoke-test on a subset of them.

    Parameters:
    mt (hl.MatrixTable): The input MatrixTable.
    is_test (bool): Flag indicating if this is a test run. Ignored if False.
    test_chrs (Optional[str]): Chromosome(s) to filter to, delimiter can be , _ or space.

    Returns:
    hl.MatrixTable: The optionally chromosome-filtered MatrixTable.
    """
    if is_test:
        test_chrs_list = re.split(r'[,_\s]+', test_chrs)
        chrs_grch38 = hl.get_reference('GRCh38').contigs
        chrs_grch38 = chrs_grch38[0:24]
        validate_input_chromosomes(test_chrs_list, chrs_grch38)
        mt = mt.filter_rows(hl.array(test_chrs_list).contains(mt.locus.contig))
        print(f"Filtering for chromosome: {test_chrs}")
    else:
        print("No chromosome filtering applied.")
    return mt


def filter_to_test_intervals_if_needed(mt: hl.MatrixTable, test_intervals: Optional[str]) -> hl.MatrixTable:
    """
    Optionally filter the MatrixTable down to one or more genomic intervals via hl.filter_intervals, e.g.
    "chr21:1-5000000" for a single small region, or a comma-delimited list of intervals (e.g. a handful
    of known hits from each of the exome/clinvar/acaf_threshold BEDs) to cheaply exercise all 3 BED types
    in one smoke test. This is the knob to use for a fast fat-box test on "a small region" rather than a
    whole chromosome -- it's a real filter on the MT/VCF data itself, unlike truncating a BED file's
    interval table, which doesn't reduce what Hail actually reads.

    Parameters:
    mt (hl.MatrixTable): The input MatrixTable.
    test_intervals (Optional[str]): A comma-delimited list of locus interval strings, or None/empty to
                                    skip filtering.

    Returns:
    hl.MatrixTable: The optionally interval-filtered MatrixTable.
    """
    if test_intervals:
        interval_strs = [s for s in re.split(r'\s*,\s*', test_intervals.strip()) if s]
        intervals = [hl.parse_locus_interval(s, reference_genome=mt.locus.dtype.reference_genome)
                    for s in interval_strs]
        mt = hl.filter_intervals(mt, intervals)
        print(f"Filtering to test interval(s): {interval_strs}")
    else:
        print("No test interval filtering applied.")
    return mt


def annotate_beds_and_filter(mt: hl.MatrixTable, bed_file_dict: Dict[str, str]) -> hl.MatrixTable:
    """
    Annotate and filter the MatrixTable based on regions specified in UCSC BED files. This script runs
    against all 3 BED files every time (exome, clinvar, acaf_threshold) -- there's one annotate/filter
    pass per bed_file_dict entry, not a choice of just one. For a test run that needs to exercise all 3
    BED types cheaply, use filter_to_test_intervals_if_needed with a few known-good intervals from each
    BED rather than subsetting the BED files themselves here.

    Parameters:
    mt (hl.MatrixTable): The input MatrixTable.
    bed_file_dict (Dict[str, str]): A dictionary mapping types of UCSC BED files to their file paths.

    Returns:
    hl.MatrixTable: The annotated and filtered MatrixTable.
    """

    def join_bed(mt: hl.MatrixTable, bed_path: str, annotation_name: str) -> hl.MatrixTable:
        """
        Parameters:
        mt (hl.MatrixTable): The Hail MatrixTable.
        bed_path (str): The gs URL to the UCSC bed file.
        annotation_name: Annotation name given to the UCSC bed file while annotating the intervals in
                         the bed file to the MatrixTable.
        Returns:
        hl.MatrixTable: An annotated MatrixTable.
        """
        ht = hl.import_bed(bed_path, reference_genome=mt.locus.dtype.reference_genome)
        print(f'UCSC Bed file at {bed_path} contains {ht.count()} intervals')
        return mt.annotate_rows(**{annotation_name: hl.is_defined(ht[mt.locus])})

    for bed_name, bed_path in bed_file_dict.items():
        mt = join_bed(mt, bed_path, bed_name)

    beds = bed_file_dict.keys()
    mt = mt.filter_rows(hl.any(*(mt[fd] for fd in beds)), keep=True)
    return mt


# this will be extracted into a shared library in the future
def make_basis_mt_from_vcf(mt: hl.MatrixTable, is_filtering_FT: bool, is_keep_as_vqsr_fields: bool,
                           max_alt_alleles: Optional[int], fields_to_drop: str) -> hl.MatrixTable:
    """
    Take a bed-annotated-and-filtered MatrixTable imported from a phased VCF and produce the same basis
    MT shape that sc1_create_basis_mt_from_vds.py produces from a VDS: recomputed variant_qc/info, with
    some fields dropped.

    Unlike the VDS path, GT and AD come directly from the VCF's own FORMAT fields -- there is no LGT/LA
    local-allele encoding to reconstruct, and no gvcf_info combiner artifact to drop. GT phase, parsed
    from the VCF's "|"-delimited genotypes, passes through untouched.

    Parameters:
    mt (hl.MatrixTable): The bed-annotated, bed-filtered input MatrixTable.
    is_filtering_FT (bool): Flag to filter GTs based on the FT field, if present. If "FT" is present and
                        does not contain "PASS", set the corresponding "GT" to missing. A real VCF's FT
                        FORMAT field, if present, is already a string (e.g. "PASS" or a semicolon-delimited
                        filter list) -- unlike the VDS/GVCF-combiner FT, which is boolean and needs to be
                        converted to a string first. That conversion step is intentionally not present here.
    is_keep_as_vqsr_fields (bool): Flag to keep VQSR fields, if an "as_vqsr" row field happens to be present.
    max_alt_alleles (Optional[int]): Maximum number of alternate alleles for filtering.
    fields_to_drop (str): Specific fields to drop, e.g., "as_vqsr, tranche_data,
                          truth_sensitivity_snp_threshold, truth_sensitivity_indel_threshold,
                          snp_vqslod_threshold, indel_vqslod_threshold", delimiter can be , _ or space.

    Returns:
    hl.MatrixTable: The resulting basis MatrixTable with recomputed variant_qc/info annotations.
    """
    if max_alt_alleles is not None:
        mt = mt.filter_rows(hl.len(mt.alleles) < max_alt_alleles)

    if is_filtering_FT and 'FT' in mt.entry:
        mt = mt.annotate_entries(GT=hl.or_missing((~hl.is_defined(mt.FT)) | (mt.FT.contains("PASS")), mt.GT))

    mt = hl.variant_qc(mt)
    computed_info = hl.struct(AC=mt.variant_qc.AC[1:], AF=mt.variant_qc.AF[1:], AN=mt.variant_qc.AN,
                              HC=mt.variant_qc.homozygote_count)
    if 'info' in mt.row:
        # A VCF-derived MT already has an `info` struct populated from the source VCF header (e.g.
        # phasing-tool annotations). Merge our recomputed fields into it instead of clobbering it, since
        # some of those original fields (e.g. phasing quality/confidence) may be needed for the eventual
        # VCF export step.
        mt = mt.annotate_rows(info=mt.info.annotate(**computed_info))
    else:
        mt = mt.annotate_rows(info=computed_info)

    mt = mt.annotate_rows(variant_qc=mt.variant_qc.drop("AC", "AF", "AN", "homozygote_count"))
    fields_to_drop_list = re.split(r'[,_\s]+', fields_to_drop)
    mt = mt.drop(*(f for f in fields_to_drop_list if f in mt.entry or f in mt.row or f in mt.col))

    if is_keep_as_vqsr_fields and ('as_vqsr' in mt.row):
        mt = mt.annotate_rows(info=mt.info.annotate(AS_VQSLOD=mt.as_vqsr.values().vqslod,
                                                     AS_YNG=mt.as_vqsr.values().yng_status))

    return mt


def create_output_path_dict(is_test: bool, test_chrs, task_label: str, version_label: str, output_gs_url: str) -> Dict[str, str]:
    """
    Generate an output path for the Hail MatrixTable based on input parameters. Mirrors
    sc1_create_basis_mt_from_vds.py so downstream SC2/SC4 stages don't need to know or care whether the
    basis MT came from a VDS or a phased VCF.

    Returns:
    Dict[str,str]: A dictionary to store the URL for the output Hail MatrixTable.
    """
    fn_append = "test_" + test_chrs if is_test else ""
    output_hail_mt_url_var = f'{output_gs_url}/{version_label}/{task_label}/{version_label}.{task_label}.{fn_append}.mt'
    print(f'Hail MT: {output_hail_mt_url_var}')
    return {"dense_mt": output_hail_mt_url_var}


def check_output_does_not_exist(path: str, overwrite_existing_outputs: bool) -> None:
    """
    Guard against silently clobbering an existing output.

    Parameters:
    path (str): The gs:// path to check.
    overwrite_existing_outputs (bool): If True, skip the check -- caller has explicitly opted in to
                                       replacing whatever is already there.

    Raises:
    FileExistsError: If `path` already exists and overwrite_existing_outputs is False.
    """
    if overwrite_existing_outputs:
        return
    if hl.hadoop_exists(path):
        raise FileExistsError(
            f"Output already exists at {path} and --overwrite_existing_outputs is false. "
            f"Pick a new --version_label/--task_label, or pass --overwrite_existing_outputs true "
            f"if you really mean to replace it."
        )


def write_dict_to_json(basis_mt_dict: Dict[str, str], basis_mt_json_file_name: str,
                       output_gs_url: str, task_label: str, version_label: str,
                       overwrite_existing_outputs: bool) -> None:
    """
    Convert a dictionary to a JSON object, write JSON object to a file, save file Google Bucket.
    """
    basis_mt_json_url = f'{output_gs_url}/{version_label}.{task_label}.{basis_mt_json_file_name}'
    check_output_does_not_exist(basis_mt_json_url, overwrite_existing_outputs)

    output_file_store_dict = "store_basis_mt_dict.json"
    with open(output_file_store_dict, "w") as output_file:
        json.dump(basis_mt_dict, output_file)
    print('finished writing json file.')
    subprocess.run(["cat", output_file_store_dict])
    no_clobber_flag = [] if overwrite_existing_outputs else ["-n"]
    subprocess.run(["gsutil", "-m", "cp"] + no_clobber_flag + [output_file_store_dict, basis_mt_json_url])
    print('finished copying json file to workspace bucket')
    print(basis_mt_json_url)
    return


# Command-line argument parsing
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process phased VCF(s) and UCSC BED files to create a filtered MatrixTable.")
    parser.add_argument('--executor_memory', help='Memory assigned to each worker', required=True)
    parser.add_argument('--executor_cores', help='CPUs assigned to each worker', required=True)
    parser.add_argument('--driver_memory', help='Memory assigned to the driver (master)', required=True)
    parser.add_argument('--reference_genome', help='Reference genome, e.g., "GRCh38"', required=True)
    parser.add_argument('--exome_bed_path', help='Path to exome UCSC BED file', required=True)
    parser.add_argument('--clinvar_bed_path', help='Path to ClinVar UCSC BED file', required=True)
    parser.add_argument('--acaf_threshold_bed_path', help='Path to aCAF threshold UCSC BED file', required=True)
    parser.add_argument('--output_gs_url', help='Base URL for the output storage location. the generated MatrixTable (MT) will be stored in the sub '
                                                'directory of this bucket.', required=True)
    parser.add_argument('--is_test', type=lambda x: (str(x).lower() == 'true'), help='Flag for test run', required=True)
    parser.add_argument('--test_chrs', help='Chromosome(s) to filter to, e.g., "chr22_chr21", only used if --vcf_url spans multiple chromosomes', required=False)
    parser.add_argument('--test_intervals', help='Comma-delimited locus interval(s) to filter to for a fast smoke '
                                                 'test, e.g. "chr21:1-5000000" or a handful of known hits from '
                                                 'each of the exome/clinvar/acaf_threshold BEDs so a test run '
                                                 'exercises all 3 BED types', required=False)
    parser.add_argument('--sample_ids', help='Comma/space-delimited list of sample IDs to subset down to, '
                                              'e.g. for a fast test on 5 samples. Omit to keep all samples.', required=False)
    parser.add_argument('--version_label', help='Version identifier for the output', required=True)
    parser.add_argument('--task_label', help='task performed in the script, e.g.: "vcf2mt"', required=True)
    parser.add_argument('--min_partitions', type=int, help='Minimum number of partitions for parallelization', required=True)
    parser.add_argument('--vcf_url', help='Path/glob, or comma/space-delimited list of paths, to the input phased VCF(s)', required=True)
    parser.add_argument('--is_filtering_FT', type=lambda x: (str(x).lower() == 'true'), help='Filter GTs based on "Failed" FT field, if present', required=True)
    parser.add_argument('--is_keep_as_vqsr_fields', type=lambda x: (str(x).lower() == 'true'), help='Keep VQSR fields, if present', required=True)
    parser.add_argument('--max_alt_alleles', type=int, help='Maximum number of alternate alleles for filtering', required=False)
    parser.add_argument('--staging_bucket', help='Temporary bucket to hold files on dataproc cluster', required=True)
    parser.add_argument('--basis_mt_json_file_name', help='File name of the basis dense MatrixTables json file, '
                                                          'suffix should be .json', required=True)
    parser.add_argument('--fields_to_drop', help='Fields to drop from the VCF-derived MT, e.g., "as_vqsr, '
                                                 'tranche_data, truth_sensitivity_snp_threshold, '
                                                 'truth_sensitivity_indel_threshold, '
                                                 'snp_vqslod_threshold, indel_vqslod_threshold", '
                                                 'delimiter can be , _ or space.', required=True)
    parser.add_argument('--overwrite_existing_outputs', type=lambda x: (str(x).lower() == 'true'), default=False,
                        help='If false (default), refuse to run if the output MT or manifest JSON already '
                             'exist. Set true to explicitly allow replacing them.', required=False)
    return parser.parse_args()


if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # All outputs (and the VCF input) must be bucket paths -- Hail writes a MatrixTable directly to a
    # bucket as it's built; there is no local-disk intermediate to copy up afterward.
    require_gs_path(args.output_gs_url, "--output_gs_url")
    for vcf_path in re.split(r'[,\s]+', args.vcf_url.strip()):
        if vcf_path:
            require_gs_path(vcf_path, "--vcf_url")

    # Initialize Hail with Spark settings
    hail_init(args.executor_memory, args.executor_cores, args.driver_memory, args.reference_genome)

    # Import the phased VCF(s) directly -- no VDS densification step needed.
    mt = import_phased_vcf(args.vcf_url, args.reference_genome, args.min_partitions)

    # Optional sample subsetting, e.g. for a fast test on a handful of samples
    mt = subset_to_samples(mt, args.sample_ids)

    # Optional chromosome filtering, for multi-chromosome test runs
    mt = perform_chr_filtering_if_needed(mt, args.is_test, args.test_chrs)

    # Optional interval filtering (hl.filter_intervals), for fast fat-box smoke tests. Pass a handful of
    # known hits from each of the 3 BEDs here to exercise all 3 BED types cheaply in one test run.
    mt = filter_to_test_intervals_if_needed(mt, args.test_intervals)

    # Create a dictionary for UCSC BED files { "exome": exome_bed_path, "clinvar": clinvar_bed_path,.. etc. }
    bed_file_dict = create_bed_dict(args.exome_bed_path, args.clinvar_bed_path, args.acaf_threshold_bed_path)

    # Annotate the MT with all 3 UCSC BED files' regions and filter it
    mt = annotate_beds_and_filter(mt, bed_file_dict)

    # Recompute variant_qc/info, apply FT/max_alt_alleles filtering, drop unneeded fields
    mt = make_basis_mt_from_vcf(mt, args.is_filtering_FT, args.is_keep_as_vqsr_fields, args.max_alt_alleles, args.fields_to_drop)

    # Generate the output path for the Hail MatrixTable
    basis_mt_dict = create_output_path_dict(args.is_test, args.test_chrs, args.task_label, args.version_label, args.output_gs_url)

    # Write basis mt dictionary to a json file as output (refuses to clobber an existing manifest
    # unless --overwrite_existing_outputs is set)
    write_dict_to_json(basis_mt_dict, args.basis_mt_json_file_name, args.output_gs_url,
                       args.task_label, args.version_label, args.overwrite_existing_outputs)

    # Export the resulting MatrixTable to the specified bucket location (refuses to clobber an
    # existing MT unless --overwrite_existing_outputs is set)
    check_output_does_not_exist(basis_mt_dict["dense_mt"], args.overwrite_existing_outputs)
    mt.write(basis_mt_dict["dense_mt"], overwrite=args.overwrite_existing_outputs)
