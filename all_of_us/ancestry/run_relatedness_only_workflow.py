"""
This script processes Variant Call Format (VCF) files along with Principal Component Analysis (PCA) projection
scores to identify related samples in genetic datasets. It creates one list:
1. A list of related samples.

NOTE:  This script does NOT generate a list of samples that can be removed from the first list using the Maximal Independent Set (MIS) method
   to account for relatedness in downstream analyses.

"""

import argparse
import logging
import hail as hl
import tarfile
import subprocess
import os


def hail_init(executor_memory: str, executor_cores: str, driver_cores: str, driver_memory: str, reference_genome: str) -> None:
    """
    Initialize Hail with specific Apache Spark configuration settings.

    Parameters:
    executor_memory (str): Memory assigned to each Spark executor (e.g., '4g').
    executor_cores (str): Number of cores assigned to each Spark executor.
    driver_cores (str): Number of cores assigned to each Spark driver.
    driver_memory (str): Memory assigned to the Spark driver node.
    reference_genome (str): Reference genome identifier (e.g., 'GRCh38').
    """
    spark_conf = {
        "spark.executor.memory": executor_memory,
        "spark.executor.cores": executor_cores,
        "spark.driver.memory": driver_memory,
        "spark.driver.cores": driver_cores

    }
    hl.init(default_reference=reference_genome, idempotent=True, spark_conf=spark_conf,
            quiet=False, skip_logging_configuration=False)


def parse_arguments() -> argparse.Namespace:
    """
    Parse and validate command-line arguments required for the script.

    Returns:
    argparse.Namespace: An object containing parsed and validated arguments.
    """
    parser = argparse.ArgumentParser(
        description="Process VCF and PCA projection files for genetic sample optimization.")
    parser.add_argument('--executor_memory', help='Memory assigned to each Spark worker')
    parser.add_argument('--executor_cores', help='CPUs assigned to each Spark worker')
    parser.add_argument('--driver_cores', help='CPUs assigned to each Spark driver')
    parser.add_argument('--driver_memory', help='Memory assigned to the Spark driver node')
    parser.add_argument('--reference_genome', help='Reference genome identifier (e.g., "GRCh38")')
    parser.add_argument('--output_gs_url', help='Output URL for the generated files')
    parser.add_argument('--task_identifier', help='Unique task identifier for the output (e.g., "aou_delta")')
    parser.add_argument('--min_partitions', type=int, help='Minimum number of partitions for Spark parallelization')
    parser.add_argument('--vcf_url', help='URL of the input VCF file')
    parser.add_argument('--pca_scores_url',
                        help='URL of the PCA scores associated with the input VCF (stored in a Hail Table)')
    parser.add_argument('--min_individual_maf', type=float, help='Minimum individual-specific minor allele frequency '
                                                                 'for relatedness')
    parser.add_argument('--statistics', help='Statistical method for relatedness (e.g., "kin")')
    parser.add_argument('--min_kinship', type=float, help='Minimum kinship threshold for identifying related samples')
    parser.add_argument('--block_size', type=int, help='Block size for matrix operations in Hail')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Initialize Hail with specified Spark settings
    hail_init(args.executor_memory, args.executor_cores, args.driver_cores, args.driver_memory, args.reference_genome)

    # Load VCF and PCA scores
    variants_mt = hl.import_vcf(args.vcf_url, force_bgz=True, min_partitions=args.min_partitions)

    # Check if the URL ends with "tar.gz"
    if args.pca_scores_url.endswith("tar.gz"):
        # Full local path for tar.gz file
        local_tar_path = "pca_scores.ht.tar.gz"
        local_ht_path = "./"

        # Download the file from GCS
        subprocess.run(['gsutil', 'cp', args.pca_scores_url, local_tar_path], check=True)

        # Decompress the file
        with tarfile.open(local_tar_path, "r:gz") as tar:
            # List all members of the tar file
            members = tar.getmembers()

            # Assuming the first member is the directory you want
            # This will get the name of the first member in the tar file
            folder_name = members[0].name.split('/')[0]
            tar.extractall()

        # Copy the decompressed file to the cluster bucket
        subprocess.run(['gsutil', '-m', 'cp', '-r', os.path.join(local_ht_path, folder_name), f'{args.output_gs_url}/'],
                       check=True)

        print(f'{args.output_gs_url}/')
        print(f'{args.output_gs_url}/{folder_name}')

        # Read the Hail table from the decompressed file in GCS
        scores_ht = hl.read_table(f'{args.output_gs_url}/{folder_name}')
    else:
        # Directly read the Hail table from the provided URL
        print(args.pca_scores_url)
        scores_ht = hl.read_table(args.pca_scores_url)

    # Identify related samples using Hail's pc_relate method
    related_samples = hl.pc_relate(variants_mt.GT, min_individual_maf=args.min_individual_maf,
                                   statistics=args.statistics,
                                   scores_expr=scores_ht[variants_mt.col_key].scores,
                                   min_kinship=args.min_kinship, block_size=args.block_size)

    related_samples = related_samples.flatten()

    related_samples.export(f"{args.output_gs_url}/{args.task_identifier}_relatedness.tsv")
    logging.info(f"Related samples generated successfully")

    logging.info(f"Relatedness outputs successfully written to: {args.output_gs_url}")
