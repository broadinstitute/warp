# Given a dense Hail MT for each bed file generated from a basis dense Hail MT, create the corresponding split
# dense Hail MatrixTables for each bed file.
# All BED files referenced in this python script are UCSC bed files (as opposed to PLINK bed files).


import argparse
import logging
import hail as hl
from typing import Dict
import json
import subprocess


# Set spark and initialize Hail
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


# Load dense mts dictionary from dense mts json file
def load_json_to_dict(input_dense_mts_json_path: str) -> Dict[str, str]:
    """
    Load a json file and output a dictionary.

    Parameters:
    input_dense_mts_json_path (str): URL of the input dense mts json file created from workflow MT2MTs.

    Returns:
    Dict[str, str]: A dictionary mapping bed dense MatrixTables to paths
    """
    subprocess.run(["gsutil", "-m", "cp", input_dense_mts_json_path, "store_dense_mts_dict.json"])
    with open("store_dense_mts_dict.json") as file:
        dense_mts_dict = json.load(file)
    print("successfully load dense_mts.json as dict.")
    return dense_mts_dict


# Get all bed names
def get_bed_names(dense_mts_dict: Dict[str, str]) -> list[str]:
    """
    Create a list of BED file names extracted from the dense MatrixTable dictionary.

    Parameters:
    dense_mts_dict (Dict[str, str]): A dictionary mapping bed dense MatrixTables to paths, generated from method "load_json_to_dict".

    Returns:
    List[str]: A list of BED names. Should match those in the VDS2MT workflow
    """
    bed_names = [k for k in sorted(list(dense_mts_dict.keys()))]
    return bed_names


# Create a dictionary mapping bed names and output split dense MTs paths
def create_output_split_dense_mts_dict(version_label: str, task_label: str,
                                       bed_names: list[str], output_gs_url: str) -> Dict[str,str]:
    """
    Generate an output Google Storage directory to store all the dense MatrixTables.

    Parameters:
    version_label (str): Version identifier for the output, example: "delta_basis_without_ext_aian_prod_gq0_"
    task_label (str): Task identifier for the output, example: "/basis_vds_exports/exports/"
    bed_names (list of str): List of the bed names, also keys of input dictionary dense_mts_dict, generated from method "get_bed_names"
    output_gs_url (str): Base URL for the output storage location, should be the workspace bucket url if working in Terra.

    Returns:
    Dict[str, str]: A dictionary mapping split bed dense MatrixTables to paths.
                    e.g., {"acaf_threshold": "gs://path-to-bucket/version/task/version.task.acaf_threshold.split.mt",
                    "clinvar": "gs://path-to-bucket/version/task/version.task.clinvar.split.mt",
                    "exome": "gs://path-to-bucket/version/task/version.task.exome.split.mt"}
    """
    split_dense_mts_dict = dict()
    for bed_name in bed_names:
        final_output_url = f'{output_gs_url}/{version_label}/{task_label}/{version_label}.{task_label}.{bed_name}.split.mt'
        split_dense_mts_dict[bed_name] = final_output_url
    return split_dense_mts_dict


# Recompute AC/AF/AN and variant_qc fields
def calculate_info(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Compute AC/AF/AN and variant_qc fields for a dense Hail MatrixTable.

    Parameters:
    mt (Hl.MatrixTable): A dense Hail MatrixTable.

    Returns:
    Hl.MatrixTable: a dense Hail MatrixTable with recomputed AC/AF/AN and variant_qc fields
    """
    fields_to_drop = ["variant_qc", "info"]
    mt = mt.drop(*(f for f in fields_to_drop if f in mt.entry or f in mt.row or f in mt.col))
    mt = mt.annotate_rows(info = hl.struct(AC=mt.variant_qc.AC[1:],
                                       AF=mt.variant_qc.AF[1:],
                                       AN=mt.variant_qc.AN,
                                       homozygote_count=mt.variant_qc.homozygote_count))
    
    # Drop duplicate nested fields that are already in INFO field rendered by call_stats()
    mt = mt.annotate_rows(variant_qc = mt.variant_qc.drop("AC", "AF", "AN", "homozygote_count"))
    return mt


# Filter variants based on field "filters", split multi_allelic variants and recompute info
def write_split_mt(mt: hl.MatrixTable, output_path: str) -> None:
    """
    Given a MatrixTable with multi-allelic variants, filter variants based on field "filters", split multi_allelic variants and recompute info.

    Parameters:
    mt (Hl.MatrixTable): A dense Hail MatrixTable.
    output_path (str): path to save the split Hail MatrixTable

    Returns:
    None
    """
    # Filter rows that have a non-missing value.  
    #  Note that this will have an issue w/ "PASS"
    if 'filters' in mt.row:
        mt = mt.filter_rows(hl.is_missing(mt.filters) | (hl.len(mt['filters']) == 0), keep=True)
    split = hl.split_multi_hts(mt)
    split = calculate_info(split)
    split.write(output_path, overwrite=True)


# Export output bed split dense MTs
def export_to_split_mt(dense_mts_dict, bed_names, split_dense_mts_dict) -> None:
    """
    Handles splitting the bed dense MatrixTables to split dense MatrixTables by bed names, and storing the MatrixTables in specified paths

    Parameters:
    dense_mts_dict (Dict[str, str]): Input dictionary mapping bed dense MatrixTables to paths.
    bed_names (list of str): List of the bed names used in the workflow VDS2MT, example ["exome", "clinvar", "acaf_threshold"]
    split_dense_mts_dict (Dict[str, str]): A dictionary mapping split bed dense MatrixTables to paths. Generated from method "create_output_split_dense_mts_dict"

    Returns:
    None
    """
    for bed_name in bed_names:
        final_output_url = split_dense_mts_dict[bed_name]
        # Input is the dense MT for the bed file.  Provided by the user (unfortunately)
        input_url = dense_mts_dict[bed_name]
        print(f'Creating Split MT for {bed_name} ====\n- Input: {input_url}\n- Output: {final_output_url}\n')
        write_split_mt(hl.read_matrix_table(input_url), final_output_url)
        print(f'Successfully wrote Split MT for {bed_name} to {final_output_url}\n')
    return


# Write dictionary to a json file as output
# This output will be used as input in the next workflow
def write_dict_to_json(split_dense_mts_dict: Dict[str, str], output_split_dense_mts_json_file_name: str,
                       output_gs_url: str, task_label: str, version_label: str) -> None:
    """
    Convert a dictionary to a JSON object and write JSON object to file.

    Parameters:
    split_dense_mts_dict: (Dict[str, str]): A dictionary mapping split bed dense MatrixTables to paths.
                    e.g., {"acaf_threshold": "gs://path-to-bucket/version/task/version.task.acaf_threshold.split.mt",
                    "clinvar": "gs://path-to-bucket/version/task/version.task.clinvar.split.mt",
                    "exome": "gs://path-to-bucket/version/task/version.task.exome.split.mt"
    output_split_dense_mts_json_file_name: (str): File name of the output json file to store the split_dense_mts_fict,
                                                suffix should be .json
    output_gs_url (str): Base URL for the output storage location, should be the workspace bucket url if working in Terra.
    task_label (str): Task identifier for the output, example: "MTs2splitMTs"
    version_label (str): Version identifier for the output, example: "delta_basis_without_ext_aian_prod_gq0_"

    Returns:
    N/A
    """
    vm_file_to_store_dict = "split_store_dense_mts_dict.json"
    with open(vm_file_to_store_dict, "w") as output_file:
        json.dump(split_dense_mts_dict, output_file)
    print('finished writing json file.')
    subprocess.run(["cat", vm_file_to_store_dict])
    split_dense_mts_json_url = f'{output_gs_url}/{version_label}/{task_label}/{version_label}.{task_label}.{output_split_dense_mts_json_file_name}'
    subprocess.run(["gsutil", "-m", "cp", vm_file_to_store_dict, split_dense_mts_json_url])
    print('finished copying json file to workspace bucket')
    return


# Command-line argument parsing
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process basis MatrixTable to create dense MatrixTable for "
                                                 "each BED name.", required=True)
    parser.add_argument('--executor_memory', help='Memory assigned to each worker', required=True)
    parser.add_argument('--executor_cores', help='CPUs assigned to each worker', required=True)
    parser.add_argument('--driver_memory', help='Memory assigned to the driver (master)', required=True)
    parser.add_argument('--reference_genome', help='Reference genome, e.g., "GRCh38"', required=True)
    parser.add_argument('--output_gs_url', help='Base URL for the output storage location. the generated MatrixTables (MTs) will be stored in the sub '
                                                'directory of this bucket.', required=True)
    parser.add_argument('--input_dense_mts_json_path', help='URL of the input dense mts json file created '
                                                      'from workflow MT2denseMTs', required=True)
    parser.add_argument('--version_label', help='Version identifier for the output, e.g.: '
                                                '"delta_basis_without_ext_aian_prod_gq0"', required=True)
    parser.add_argument('--task_label', help='taks performed in the script, e.g.: denseMTs2splitMTs"', required=True)
    parser.add_argument('--staging_bucket', help='Temporary bucket to hold files on dataproc cluster', required=True)
    parser.add_argument('--output_split_dense_mts_json_file_name', help='File name of the bed split dense MatrixTables json file, '
                                                                 'suffix should be .json', required=True)
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
    dense_mts_dict = load_json_to_dict(args.input_dense_mts_json_path)

    # Get list of bed names
    bed_names = get_bed_names(dense_mts_dict)

    # Create a dictionary mapping bed names and output split dense MTs paths
    split_dense_mts_dict = create_output_split_dense_mts_dict(args.version_label, args.task_label,
                                                              bed_names, args.output_gs_url)

    # Export bed split MTs to bucket
    export_to_split_mt(dense_mts_dict, bed_names, split_dense_mts_dict)

    # Write dense mts dictionary to a json file as output
    write_dict_to_json(split_dense_mts_dict, args.output_split_dense_mts_json_file_name, args.output_gs_url,
                       args.task_label, args.version_label)
