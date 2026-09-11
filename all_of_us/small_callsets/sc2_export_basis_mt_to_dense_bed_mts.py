# Given a basis Hail MatrixTable created from VDS and annotated by UCSC BED names using
# the workflow "sc1_create_basis_mt_from_vds" (https://github.com/broadinstitute/aou_small_callsets/blob/main/sc1_create_basis_mt_from_vds.wdl),
# create corresponding dense Hail MatrixTables for each UCSC BED file.
# All UCSC BED files referenced in this python script are UCSC BED files (as opposed to PLINK bed files).

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


# Load basis dense mt dictionary from basis dense json file
def load_json_to_dict(input_basis_dense_mt_json_path: str) -> Dict[str, str]:
    """
    Load a json file and output a dictionary.

    Parameters:
    input_basis_dense_mt_json_path (str): URL of the input basis dense mt json file created from workflow VDS2MT.
                                    The key should be "dense_mt".

    Returns:
    Dict[str, str]: A dictionary mapping UCSC BED dense MatrixTables to paths.
    """
    json_file_name_in_vm = "store_basis_dense_mt_dict.json"
    subprocess.run(["gsutil", "-m", "cp", input_basis_dense_mt_json_path, json_file_name_in_vm])
    with open(json_file_name_in_vm) as file:
        basis_dense_mt_dict = json.load(file)
    print("successfully load basis dense json as dict.")
    return basis_dense_mt_dict


# Get UCSC BED names according to their corresponding paths
# Still use UCSC BED files paths as input to make sure the UCSC BED names are correct
def get_bed_names(exome_bed_path: str, clinvar_bed_path: str, acaf_threshold_bed_path: str) -> list[str]:
    """
    Create a list of UCSC BED file names according to their file paths.
    
    Parameters:
    exome_bed_path (str): Path to the exome UCSC BED file.
    clinvar_bed_path (str): Path to the ClinVar UCSC BED file.
    acaf_threshold_bed_path (str): Path to the ACAF threshold UCSC BED file.
    
    Returns:
    list[str]: A list of UCSC BED files according to their corresponding paths.
    """
    bed_file_dict = {"exome": exome_bed_path,
                    "clinvar": clinvar_bed_path,
                    "acaf_threshold": acaf_threshold_bed_path}
    bed_names = [k for k in sorted(list(bed_file_dict.keys()))]
    return bed_names


# Filter a mt to a named UCSC BED interval filter
def filter_dense_to_bed(mt: hl.MatrixTable, name: str, bed_names: list[str]) -> hl.MatrixTable:
    """
    Handles filtering the basis MatrixTable to a dense MatrixTable based on the UCSC BED name.
    The MatrixTable includes flags indicating which bed files contain a given variant by adding an
    attribute to each variant with the bed file key (e.g., 'exome').
    
    Parameters:
    mt: (hl.MatrixTable): The input basis MatrixTable created from VDS, annotated by UCSC BED names
    name: (str): name of the UCSC BED file: exome or clinvar or acaf_threshold
    bed_names (list of str): List of the UCSC BED names used in the workflow "sc1_create_basis_mt_from_vds"
                            https://github.com/broadinstitute/aou_small_callsets/blob/main/sc1_create_basis_mt_from_vds.wdl ,
                            this argument is the output of method "get_bed_names".
    """

    return mt.filter_rows(mt[name]).drop(*bed_names)


# Create a dictionary mapping UCSC BED names with output dense mts paths
# No need to add test component in this script as if_test is determined
# in the first workflow sc1_create_basis_mt_from_vds
def create_output_dense_mts_dict(version_label: str, task_label: str, bed_names: list[str],
                                 output_gs_url: str) -> Dict[str,str]:
    """
    Generate an output Google Storage directory to store all the dense MatrixTables.
    
    Parameters:
    version_label (str): Version identifier for the output, example: "delta_basis_without_ext_aian_prod_gq0_"
    task_label (str): Task identifier for the output, example: "/basis_vds_exports/exports/"
    bed_names (list of str): List of the UCSC BED names used in the workflow "sc1_create_basis_mt_from_vds"
                            https://github.com/broadinstitute/aou_small_callsets/blob/main/sc1_create_basis_mt_from_vds.wdl ,
                            this argument is the output of method "get_bed_names".
    output_gs_url (str): Base URL for the output storage location, should be the workspace bucket url if working in Terra, must have trailing in the end
    
    Returns:
    Dict[str, str]: A dictionary mapping UCSC BED dense MatrixTables to paths.
    """
    identifier = f'{version_label}{len(bed_names)}regions'
    # MUST have trailing '/'
    output_bucket = f'{output_gs_url}{version_label}{task_label}/'
    out_dense_mt_base_url = output_bucket + identifier

    dense_mts_dict = {}
    for bed_name in bed_names:
        output_url = f'{output_gs_url}/{version_label}.{task_label}.{bed_name}.mt'
        dense_mts_dict[bed_name] = output_url
        print(f'adding {bed_name} mt path to dict.')
    return dense_mts_dict
    


# Filter basis MT and write to dense mts by UCSC BED names
def filter_mt_by_bed_name(bed_names: list, dense_mts_dict: Dict[str,str], basis_dense_mt_dict: Dict) -> None:
    """
    Handles filtering the basis MatrixTable to dense MatrixTables by UCSC BED names, and storing the MatrixTables in specified paths
    
    Parameters:
    dense_mts_dict (Dict[str, str]): A dictionary mapping UCSC BED dense MatrixTables to paths, generated from method "create_output_dense_mt_base_url".
                                    The MatrixTable should have a row field `bed_names`, with values from the list "bed_names".
    bed_names (list of str): List of the UCSC BED names used in the workflow "sc1_create_basis_mt_from_vds"
                            https://github.com/broadinstitute/aou_small_callsets/blob/main/sc1_create_basis_mt_from_vds.wdl ,
                            this argument is the output of method "get_bed_names".
    basis_dense_mt_dict (dict): Dict of basis dense mt name and path created from vds, the key should be "dense_mt".

    Returns:
    N/A
    """
    for bed_name in bed_names:
        output_url = dense_mts_dict[bed_name]
        logging.info(f'Starting {bed_name.upper()} dense MT (multi): {output_url}')
        print(f'Starting {bed_name.upper()} dense MT (multi): {output_url}')
        mt_url = basis_dense_mt_dict["dense_mt"]
        dense_mt_bed = filter_dense_to_bed(hl.read_matrix_table(mt_url), bed_name)
        dense_mt_bed.write(output_url, overwrite = True)
    return


# Write dictionary to a json file as output
# This output will be used as input in the next workflow
def write_dict_to_json(dense_mts_dict: Dict[str, str], output_dense_mts_json_file_name: str,
                       output_gs_url: str, task_label: str, version_label: str) -> None:
    """
    Convert dictionary to JSON object and write JSON object to file.

    Parameters:
    dense_mts_dict: (Dict[str, str]): A dictionary mapping UCSC BED dense MatrixTables to paths, generated from method "create_output_dense_mt_base_url"
    output_dense_mts_json_file_name: (str): File name of the UCSC BED dense MatrixTables json file, suffix should be .json

    Returns:
    N/A
    """
    dens_mts_json_file_vm_name = "store_dense_mts_dict.json"
    with open(dens_mts_json_file_vm_name, "w") as output_file:
        json.dump(dense_mts_dict, output_file)
    print('finished writing json file.')
    subprocess.run(["cat", dens_mts_json_file_vm_name])
    dense_mts_json_url = f'{output_gs_url}/{version_label}.{task_label}.{output_dense_mts_json_file_name}'
    subprocess.run(["gsutil", "-m", "cp", dens_mts_json_file_vm_name, dense_mts_json_url])
    print('finished copying json file to workspace bucket')
    return

    
# Command-line argument parsing
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process basis MatrixTable to create dense MatrixTable for each UCSC BED name.")
    parser.add_argument('--executor_memory', help='Memory assigned to each worker', required=True)
    parser.add_argument('--executor_cores', help='CPUs assigned to each worker', required=True)
    parser.add_argument('--driver_memory', help='Memory assigned to the driver (master)', required=True)
    parser.add_argument('--reference_genome', help='Reference genome, e.g., "GRCh38"', required=True)
    parser.add_argument('--exome_bed_path', help='Path to exome UCSC BED file', required=True)
    parser.add_argument('--clinvar_bed_path', help='Path to ClinVar UCSC BED file', required=True)
    parser.add_argument('--acaf_threshold_bed_path', help='Path to aCAF threshold UCSC BED file', required=True)
    parser.add_argument('--output_gs_url', help='Base URL for the output storage location. the generated MatrixTable (MT) will be stored in the sub '
                                                'directory of this bucket.', required=True)
    parser.add_argument('--version_label', help='Version identifier for the output, e.g.: '
                                                '"delta_basis_without_ext_aian_prod_gq0_"', required=True)
    parser.add_argument('--task_label', help='task performed in the script, e.g.: "mt2denseMTs"', required=True)
    parser.add_argument('--staging_bucket', help='Temporary bucket to hold files on dataproc cluster', required=True)
    parser.add_argument('--output_dense_mts_json_file_name', help='File name of the UCSC BED dense MatrixTables json file, '
                                                           'suffix should be .json', required=True)
    parser.add_argument('--input_basis_dense_mt_json_path', help='URL of the input basis dense mt json file created from '
                                                                 'workflow sc1_filter_vds_to_basis_MT_by_bed. The key '
                                                                 'should be "dense_mt"', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Initialize Hail with Spark settings
    # The Hail team will likely make this the default behavior in a later version.
    hail_init(args.executor_memory, args.executor_cores, args.driver_memory, args.reference_genome)
    
    # Use the new shuffle that handles preemptible nodes better.
    hl._set_flags(use_new_shuffle='1')
    
    # Get list of UCSC BED names
    bed_names = get_bed_names(args.exome_bed_path, args.clinvar_bed_path, args.acaf_threshold_bed_path)
    
    # Create a dictionary mapping output bed dense MatrixTables to paths
    dense_mts_dict = create_output_dense_mts_dict(args.version_label, args.task_label,
                                                  bed_names, args.output_gs_url)

    # Load basis dense mt dictionary from basis dense json file
    basis_dense_mt_dict = load_json_to_dict(args.input_basis_dense_mt_json_path)

    # Filter basis MT and export to dense mts by UCSC BED names
    filter_mt_by_bed_name(bed_names, dense_mts_dict, basis_dense_mt_dict)


    # Write dense mts dictionary to a json file as output
    write_dict_to_json(dense_mts_dict, args.output_dense_mts_json_file_name, args.output_gs_url,
                       args.task_label, args.version_label)
