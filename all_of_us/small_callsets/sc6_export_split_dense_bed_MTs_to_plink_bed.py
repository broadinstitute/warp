# Given a split dense MatrixTable with multi-allelic variants for each bed name,
# create corresponding PLINK BED files by chromosome for each UCSC BED file.

# PLINK BED is not the same as UCSC BED. UCSC BED were used to generate the input intervals. PLINK bed is a variant file format.

# PLINK BED outputs are found in the specified gs bucket {output_gs_url}/{version_label}/{task_label_label}/{bed_name}/{version_label}.{task_label_label}.{bed_name}.{c}


#export split dense MatrixTable to plink bed files
import hail as hl
import argparse
import pandas as pd
from typing import Dict
import subprocess
import json
import re


# Set spark and initialize Hail
def hail_init(executor_memory: str, executor_cores: str, driver_memory: str, reference_genome: str) -> None:
    """
    Initialize Hail with specific Spark configuration settings.

    Parameters:
    executor_memory (str): Memory assigned to each Spark executor, e.g., "15g".
    executor_cores (str): Number of cores assigned to each Spark executor, e.g., "8".
    driver_memory (str): Memory assigned to the Spark driver, e.g., "200".
    reference_genome (str): Reference genome version_label, e.g., "GRCh38".
    """
    spark_conf = {
        "spark.executor.memory": executor_memory,
        "spark.executor.cores": executor_cores,
        "spark.driver.memory": driver_memory
    }
    hl.init(default_reference=reference_genome, idempotent=True, spark_conf=spark_conf,
            quiet=True, skip_logging_configuration=True)


# Load split dense mts dictionary from dense mts json file
def load_json_to_dict(input_split_dense_mts_json_path: str) -> Dict[str, str]:
    """
    Load a json file and output a dictionary.

    Parameters:
    input_split_dense_mts_json_path (str): URL of the input dense mts json file.

    Returns:
    Dict[str, str]: A dictionary mapping bed dense MatrixTables to paths
    """
    subprocess.run(["gsutil", "-m", "cp", input_split_dense_mts_json_path, "split_store_dense_mts_dict.json"])
    with open("split_store_dense_mts_dict.json") as file:
        split_dense_mts_dict = json.load(file)
    print("successfully load dense_mts.json as dict.")
    return split_dense_mts_dict


# Get bed names
def get_bed_names(split_dense_mts_dict: Dict[str, str]) -> list[str]:
    """
    Create a list of BED file names extracting from the dense MatrixTable dictionary.

    Parameters:
    split_dense_mts_dict (Dict[str, str]): A dictionary mapping bed dense MatrixTables to paths, generated from method "load_json_to_dict".

    Returns:
    List[str]: A list of all BED names
    """
    bed_names = [k for k in sorted(list(split_dense_mts_dict.keys()))]
    return bed_names


# Export dense Hail MatrixTable to plink bed files
# Requires split mt as input
# write split dense MT to plink bed files by chromosome
def write_chr_plink(mt: hl.MatrixTable, chrom: str, base_path: str) -> None:
    """
    Export a Hail MatrixTable to PLINK BED files by chromosome.

    Parameters:
    mt (hl.MatrixTable): A dense Hail MatrixTable.
    chrom (str): chromosome in GRCh38, e.g., 'chr22'.
    base_path (str): base path of output plink files.

    Returns:
    None
    """
    print(f'writing PLINK for chr={chrom}')
    mt = mt.filter_rows(mt.locus.contig == chrom)
    hl.export_plink(mt, base_path, ind_id = mt.s)
    return


# Write a list to file
def write_output_list(output_list: list[str], out_gs_url: str) -> None:
    """
    Write a list of output file paths to a file in google bucket, with one line per file and no headers.

    Parameters:
    output_list (list): A list of paths to outputs.
    out_gs_url (str): path to the file to store the list.

    Returns:
    None
    """
    df = pd.Series(output_list)
    df.to_csv(out_gs_url, sep="\t", index=None, header=None)
    return


# Determine bucket to store outputs
def determine_export_bucket(output_bucket: str, bed_name: str, export_type: str) -> str:
    """
    Determine the Google bucket to store the outputs.

    Parameters:
    output_bucket (str): URL of the bucket.
    bed_name (str): One of the three bed names: 'acaf_threshold', 'clinvar', 'exome'.
    export_type (str): The file type to export. e.g., 'plink_bed'

    Returns:
    Str: the base URL of the Google bucket to store the outputs.
    """
    output_export_bucket = f'{output_bucket}{bed_name}/{export_type}/'
    return output_export_bucket


# Get bed names to skip analyzing
def skip_beds(skip_beds_names) -> list[str]:
    """
    Get the list of bed names to analyze.

    Parameters:
    skip_beds_names (str): String to specify the bed names, e.g. 'exome,chinvar', delimiter can be',' or '_', or space.

    Returns:
    List: A list of bed names, e.g. ['exome', 'clinvar']
    """
    if skip_beds_names == "None":
        skip_beds_list=[]
    else:
        skip_beds_list = re.split(r'[,_\s]+', skip_beds_names)
    return skip_beds_list


# Upload output list to manifest file
def upload_manifests(final_outputs: str, extension: str, output_export_bucket: str) -> None:
    """
    Write the URL of outputs to a manifest file in Google Bucket.

    Parameters:
    final_outputs (str): URL of the output file.
    extension (str): Type of the file, e.g., plink_bed.
    output_export_bucket (str): Base URL to store the manifest files.

    Returns:
    None
    """
    output_list = "\n".join(final_outputs)
    print(output_list)
    write_output_list(final_outputs, f'{output_export_bucket}manifest.{extension}.txt')
    print("------------")
    return


# Get the list of chromosomes to work in the analysis
def get_chrs_list(chrs: str) -> list[str]:
    """
    Get the list of chromosomes to analyze.

    Parameters:
    chrs (str): Chromosomes to analyze, e.g., 'chr22,chr21' or 'chr22_chr21', delimiter can be , _ or space.

    Returns:
    List: A list of chromosomes to analyze.
    """
    if chrs == "None":
        chrs_list = hl.get_reference('GRCh38').contigs
    # We only support the first 24 contigs (chr1-22, x, y)
        chrs_list = chrs_list[0:24]
        print(chrs_list)
    else:
        chrs_list = re.split(r'[,_\s]+', chrs)
        print(chrs_list)
    return chrs_list


# Write to plink bed by chr
def export_to_plinkbed(output_gs_url: str, bed_names: str, split_dense_mts_dict: Dict[str,str],
                       skip_beds_list: list, chrs_list: list, version_label: str, task_label: str) -> None:
    """
    Export a dense Hail MatrixTable to Plink BED files.

    Parameters:
    output_gs_url: Basis URL of the Google bucket for outputs.
    bed_names: A list of BED names. Should match those in the VDS2MT workflow.
    split_dense_mts_dict: (Dict[str, str]): A dictionary mapping bed dense MatrixTables to paths, generated from method "load_json_to_dict".
    skip_beds_list: (list): A list of bed names, e.g. ['exome', 'clinvar'].
    chrs_list (list): A list of chromosomes to analyze.
    version_label (str): Version of the whole analysis, e.g.: "delta_basis_without_ext_aian_prod_gq0_3regions"')
    task_label (str): Task performed in the script, e.g.: "mt2plink_bed"

    Returns:
    None
    """
    output_bucket = f'{output_gs_url}/{version_label}/{task_label}/'
    # Leave this list empty if you want to process all split MTs.  Otherwise, list ones to skip.
    #  This can be used to resume failures, but only at the granularity of a corpus of VCF files.
    for bed_name in bed_names:
        if bed_name in skip_beds_list:
            print(f'Told to skip {bed_name}...')
            continue
        final_outputs = []
        output_export_bucket = determine_export_bucket(output_bucket, bed_name, "plink_bed")
        print("Export bucket: " + output_export_bucket)

        for c in chrs_list:
            print(f'Starting {bed_name.upper()} plink bed: {c}')

            # This is the basename.  Output will be three files:  `{final_output}.bed`, `{final_output}.bim`, and `{final_output}.fam`
            final_output = f'{output_export_bucket}{version_label}.{bed_name}.{c}'

            split_mt_url = split_dense_mts_dict[bed_name]
            write_chr_plink(hl.read_matrix_table(split_mt_url), c, final_output)

            print(f'{c}:  {final_output}.bed')
            print(f'{c}:  {final_output}.bim')
            print(f'{c}:  {final_output}.fam')
            final_outputs.append([final_output + ".bed", final_output + ".bim", final_output + ".fam"])


        upload_manifests([f[0] for f in final_outputs], "bed", output_export_bucket)
        upload_manifests([f[1] for f in final_outputs], "bim", output_export_bucket)
        upload_manifests([f[2] for f in final_outputs], "fam", output_export_bucket)
    return

# Command-line argument parsing
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process basis MatrixTable to create dense MatrixTable "
                                                 "for each BED name.", required=True)
    parser.add_argument('--executor_memory', help='Memory assigned to each worker', required=True)
    parser.add_argument('--executor_cores', help='CPUs assigned to each worker', required=True)
    parser.add_argument('--driver_memory', help='Memory assigned to the driver (master)', required=True)
    parser.add_argument('--reference_genome', help='Reference genome, e.g., "GRCh38"', required=True)
    parser.add_argument('--output_gs_url', help='Base URL for the output storage location. the generated plink bed files will be stored in the sub '
                                                'directory of this bucket.', required=True)
    parser.add_argument('--input_split_dense_mts_json_path', help='URL of the input split dense mts json file '
                                                            'created from workflow denseMTs2splitMTs', required=True)
    parser.add_argument('--version_label', help='Version version_label for the output, e.g.: '
                                                '"delta_basis_without_ext_aian_prod_gq0_3regions"', required=True)
    parser.add_argument('--task_label', help='task_label performed in the script, e.g.: '
                                             '"splitMTs2plinkBED"', required=True)
    parser.add_argument('--skip_beds_names', help='Bed names to skip if not working on all beds,e.g., '
                                                  '"exome,clinvar", or "None", no space after comma', required=True)
    parser.add_argument('--chrs', help='Chromosomes to work on. e.g., "chr22,chr21", '
                                       'no space after comma', required=True)
    parser.add_argument('--staging_bucket', help='Temporary bucket to hold files on dataproc cluster', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Initialize Hail with Spark settings
    # The Hail team will likely make this the default behavior in a later version.
    hail_init(args.executor_memory, args.executor_cores, args.driver_memory, args.reference_genome)

    # Use the new shuffle that handles preemptible nodes better.
    hl._set_flags(use_new_shuffle='1')

    # Get dense mts dict
    split_dense_mts_dict = load_json_to_dict(args.input_split_dense_mts_json_path)

    # Get list of bed names
    bed_names = get_bed_names(split_dense_mts_dict)

    # Get skip beds list
    skip_beds_list = skip_beds(args.skip_beds_names)

    #Get list of chromosomes to work on
    chrs_list = get_chrs_list(args.chrs)

    # Write to plink bed by chromosome
    export_to_plinkbed(args.output_gs_url, bed_names, split_dense_mts_dict, skip_beds_list,
                       chrs_list, args.version_label, args.task_label)
