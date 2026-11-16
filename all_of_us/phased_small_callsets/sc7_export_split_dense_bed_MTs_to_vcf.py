# Given a split dense Hail MT (post sc4_export_dense_bed_MTs_to_split_dense_MTs.py), export it as a VCF
# per callset/chromosome. This is the first time this pipeline produces a VCF output -- there is no
# equivalent VDS-input script to mirror.
#
# Reads SC4's split-dense-MT json manifest and writes one VCF per BED name. Output granularity is
# per-chromosome by way of the pipeline running per-chromosome end to end (one chromosome's phased VCF in
# at SC1, one chromosome's data throughout) -- this script does not split by contig itself.
#
# STATUS: not yet run. Open items:
#   - SC4's calculate_info() rebuilds `info` from scratch as {AC, AF, AN, homozygote_count}, dropping
#     GQ/AVSAD/HC. build_phased_vcf_metadata() below still declares GQ/AVSAD/HC, so the declared header
#     and the actual data have diverged. Pending the GQ/AVSAD/HC-vs-homozygote_count decision, this
#     metadata may need trimming, or SC4 may need to stop clobbering `info`.
#   - Whether the LowQual/ExcessHet/NO_HQ_GENOTYPES FILTER definitions are ever actually needed
#     (SC4 drops any row with a non-empty `filters` field before splitting)
#   - Whether `##reference=` is needed by anything downstream (confirmed Hail does NOT emit one, and
#     offers no metadata slot for it -- would require post-processing)
#   - Indexing is deliberately NOT done here. Per guidance, index separately with tabix/bcftools as its
#     own cheap non-Hail task rather than Hail's tabix=True.
#
# Deviation from the sibling scripts worth knowing: this one does not accept --staging_bucket. Every
# sibling parses it and then never references it; rather than carry a vestigial required argument, it is
# omitted here and correspondingly not passed by the WDL.
#
# See all_of_us/phasing/AoU_VCF_v9.py and all_of_us/phasing/remove_phased_samples_X.py
# (build_default_metadata) for the reference conventions this is based on.


import argparse
import hail as hl
from typing import Dict
import json
import subprocess


def hail_init(executor_memory: str, executor_cores: str, driver_memory: str, reference_genome: str) -> None:
    """
    Initialize Hail with specific Spark configuration settings.

    Parameters:
    executor_memory (str): Memory assigned to each Spark executor, e.g., "10g".
    executor_cores (str): Number of cores assigned to each Spark executor, e.g., "4".
    driver_memory (str): Memory assigned to the Spark driver, e.g., "80g".
    reference_genome (str): Reference genome identifier, e.g., "GRCh38".
    """
    spark_conf = {
        "spark.executor.memory": executor_memory,
        "spark.executor.cores": executor_cores,
        "spark.driver.memory": driver_memory
    }
    hl.init(default_reference=reference_genome, idempotent=True, spark_conf=spark_conf,
            quiet=True, skip_logging_configuration=True)


def load_json_to_dict(input_split_dense_mts_json_path: str) -> Dict[str, str]:
    """
    Load SC4's split dense MTs manifest and return it as a dictionary.

    Parameters:
    input_split_dense_mts_json_path (str): gs:// URL of the json file written by SC4.

    Returns:
    Dict[str, str]: A dictionary mapping bed name to split dense MatrixTable path.
    """
    json_file_name_in_vm = "store_split_dense_mts_dict.json"
    subprocess.run(["gsutil", "-m", "cp", input_split_dense_mts_json_path, json_file_name_in_vm])
    with open(json_file_name_in_vm) as file:
        split_dense_mts_dict = json.load(file)
    print("successfully loaded split_dense_mts json as dict.")
    return split_dense_mts_dict


def get_bed_names(split_dense_mts_dict: Dict[str, str]) -> list[str]:
    """
    Create a sorted list of BED names from the manifest's keys.

    Parameters:
    split_dense_mts_dict (Dict[str, str]): Manifest dict from load_json_to_dict.

    Returns:
    list[str]: Sorted BED names, e.g. ["acaf_threshold", "clinvar", "exome"].
    """
    return [k for k in sorted(list(split_dense_mts_dict.keys()))]


def create_output_vcfs_dict(version_label: str, task_label: str, bed_names: list[str],
                            output_gs_url: str) -> Dict[str, str]:
    """
    Build the output VCF path for each BED name.

    Returns:
    Dict[str, str]: A dictionary mapping bed name to output .vcf.bgz path. The .bgz suffix matters --
                    hl.export_vcf block-gzips based on the file extension.
    """
    vcfs_dict = dict()
    for bed_name in bed_names:
        vcfs_dict[bed_name] = (f'{output_gs_url}/{version_label}/{task_label}/'
                               f'{version_label}.{task_label}.{bed_name}.vcf.bgz')
    return vcfs_dict


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


def export_all_to_vcf(split_dense_mts_dict: Dict[str, str], bed_names: list[str],
                      vcfs_dict: Dict[str, str], overwrite_existing_outputs: bool) -> None:
    """
    Export each BED name's split dense MatrixTable to its own VCF.

    Parameters:
    split_dense_mts_dict (Dict[str, str]): bed name -> input split dense MT path (from SC4's manifest).
    bed_names (list[str]): BED names to process, from get_bed_names.
    vcfs_dict (Dict[str, str]): bed name -> output VCF path, from create_output_vcfs_dict.
    overwrite_existing_outputs (bool): If False, refuse to overwrite an existing VCF.

    Returns:
    None
    """
    for bed_name in bed_names:
        input_url = split_dense_mts_dict[bed_name]
        output_url = vcfs_dict[bed_name]
        print(f'Exporting VCF for {bed_name} ====\n- Input: {input_url}\n- Output: {output_url}\n')
        export_mt_to_vcf(hl.read_matrix_table(input_url), output_url, overwrite_existing_outputs)
        print(f'Successfully wrote VCF for {bed_name} to {output_url}\n')
    return


def write_dict_to_json(vcfs_dict: Dict[str, str], output_vcfs_json_file_name: str,
                       output_gs_url: str, task_label: str, version_label: str,
                       overwrite_existing_outputs: bool) -> None:
    """
    Write the bed-name -> VCF path mapping to a json manifest in the output bucket, for the downstream
    indexing step to consume.

    Returns:
    None
    """
    vcfs_json_url = (f'{output_gs_url}/{version_label}/{task_label}/'
                     f'{version_label}.{task_label}.{output_vcfs_json_file_name}')
    check_output_does_not_exist(vcfs_json_url, overwrite_existing_outputs)

    local_file = "store_vcfs_dict.json"
    with open(local_file, "w") as output_file:
        json.dump(vcfs_dict, output_file)
    print('finished writing json file.')
    subprocess.run(["cat", local_file])
    no_clobber_flag = [] if overwrite_existing_outputs else ["-n"]
    subprocess.run(["gsutil", "-m", "cp"] + no_clobber_flag + [local_file, vcfs_json_url])
    print('finished copying json file to workspace bucket')
    print(vcfs_json_url)
    return


# Command-line argument parsing
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Export split dense MatrixTables to per-BED VCFs.")
    parser.add_argument('--executor_memory', help='Memory assigned to each worker', required=True)
    parser.add_argument('--executor_cores', help='CPUs assigned to each worker', required=True)
    parser.add_argument('--driver_memory', help='Memory assigned to the driver (master)', required=True)
    parser.add_argument('--reference_genome', help='Reference genome, e.g., "GRCh38"', required=True)
    parser.add_argument('--input_split_dense_mts_json_path', help='gs:// URL of the split dense MTs json '
                                                                  'manifest written by SC4', required=True)
    parser.add_argument('--output_gs_url', help='Base gs:// URL for the output VCFs', required=True)
    parser.add_argument('--version_label', help='Version identifier for the output', required=True)
    parser.add_argument('--task_label', help='Task performed in the script, e.g.: "splitmts2vcf"', required=True)
    parser.add_argument('--output_vcfs_json_file_name', help='File name of the output VCFs json manifest, '
                                                             'suffix should be .json', required=True)
    parser.add_argument('--overwrite_existing_outputs', type=lambda x: (str(x).lower() == 'true'), default=False,
                        help='If false (default), refuse to run if an output VCF or the manifest json '
                             'already exists. Set true to explicitly allow replacing them.', required=False)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Outputs must be bucket paths -- Hail writes directly to GCS, there is no local intermediate.
    require_gs_path(args.output_gs_url, "--output_gs_url")
    require_gs_path(args.input_split_dense_mts_json_path, "--input_split_dense_mts_json_path")

    hail_init(args.executor_memory, args.executor_cores, args.driver_memory, args.reference_genome)

    # Read SC4's manifest to find the per-BED split dense MTs
    split_dense_mts_dict = load_json_to_dict(args.input_split_dense_mts_json_path)
    bed_names = get_bed_names(split_dense_mts_dict)

    # Build output VCF paths and export one VCF per BED name
    vcfs_dict = create_output_vcfs_dict(args.version_label, args.task_label, bed_names, args.output_gs_url)
    export_all_to_vcf(split_dense_mts_dict, bed_names, vcfs_dict, args.overwrite_existing_outputs)

    # Write the manifest of VCF paths for the downstream indexing step
    write_dict_to_json(vcfs_dict, args.output_vcfs_json_file_name, args.output_gs_url,
                       args.task_label, args.version_label, args.overwrite_existing_outputs)
