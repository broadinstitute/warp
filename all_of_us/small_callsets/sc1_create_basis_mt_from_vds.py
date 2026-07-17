# Given a VDS and a list of bed file, create a matrixTable of filtered variants to regions in these bed files.
# All bed files referenced in this python script are UCSC bed files (as opposed to PLINK bed files).


import argparse
import hail as hl
from typing import Dict, Optional
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


def create_output_path_dict(is_test: bool, test_chrs, task_label: str, version_label: str, output_gs_url: str) -> Dict[str,str]:
    """
    Generate an output path for the Hail MatrixTable based on input parameters.

    Parameters:
    is_test (bool): Flag indicating if this is a test run. A test run will only work on chromosomes specified in the
                    argument "test_chrs" in the workflow. It will work on all chr1-Y chromosomes if it is not a test run.
    test_chrs (str): Specific chromosomes to filter, if applicable, for example "chr21,chr22", delimiter can be , _ or space.
    version_label (str): Version identifier for the output. It is the same identifier for the
                        callset and remains the same for all workflows in this repository.
    bed_file_dict (Dict[str, str]): Dictionary of UCSC BED file types and paths.
    output_gs_url (str): Base URL for the output storage location. All files will be stored in the sub directories
                        in this bucket.

    Returns:
    Dict[str,str]: A dictionary to store the URL for the output Hail MatrixTable.
    """
    fn_append = "test_" + test_chrs if is_test else ""
    output_hail_mt_url_var = f'{output_gs_url}/{version_label}/{task_label}/{version_label}.{task_label}.{fn_append}.mt'
    print(f'Hail MT: {output_hail_mt_url_var}')
    basis_mt_dict = {}
    basis_mt_dict["dense_mt"] = output_hail_mt_url_var
    return basis_mt_dict

# Write dictionary to a json file as output
# This output will be used as input in the next workflow
def write_dict_to_json(basis_mt_dict: Dict[str, str], basis_mt_json_file_name: str,
                       output_gs_url: str, task_label: str, version_label: str) -> None:
    """
    Convert a dictionary to a JSON object, write JSON object to a file, save file Google Bucket.

    Parameters:
    basis_mt_dict (Dict[str, str]): A dictionary mapping UCSC BED dense MatrixTables to paths, generated from method "create_output_dense_mt_base_url".
    basis_mt_json_file_name (str): File name of the UCSC BED dense MatrixTables json file, suffix should be .json.
    output_gs_url (str): Base URL for the output storage location. All files will be stored in the sub directories
                        in this bucket.
    task_label (str): Task identifier for the output, example: "vds2mt".
    version_label(str): Version identifier for the output. It is the same identifier for the
                        callset and remains the same for all workflows in this repository.
    Returns:
    None
    """
    output_file_store_dict = "store_basis_mt_dict.json"
    with open(output_file_store_dict, "w") as output_file:
        json.dump(basis_mt_dict, output_file)
    print('finished writing json file.')
    subprocess.run(["cat", output_file_store_dict])
    basis_mt_json_url = f'{output_gs_url}/{version_label}.{task_label}.{basis_mt_json_file_name}'
    subprocess.run(["gsutil", "-m", "cp", output_file_store_dict, basis_mt_json_url])
    print('finished copying json file to workspace bucket')
    print(basis_mt_json_url)
    return


def perform_repartition_if_needed(min_partitions: int, vds_url: str) -> hl.vds.VariantDataset:
    """
    Handle repartitioning of the Variant Dataset (VDS) for optimization purposes.

    Parameters:
    min_partitions (int): Minimum number of partitions for optimization.
    vds_url (str): URL of the input VDS.

    Returns:
    hl.vds.VariantDataset: The optionally repartitioned VDS.
    """
    print(vds_url)
    vds = hl.vds.read_vds(vds_url)
    vds_ref_partitions = vds.reference_data.n_partitions()
    vds_var_partitions = vds.variant_data.n_partitions()
    print(f"Reference partitions: {vds_ref_partitions}")
    print(f"Variant partitions: {vds_var_partitions}")

    is_need_repartition = (vds_ref_partitions != vds_var_partitions) or \
                       (vds_ref_partitions < min_partitions) or \
                       (vds_var_partitions < min_partitions)

    if is_need_repartition:
        print("Repartitioning is required.")
        intervals = vds.reference_data._calculate_new_partitions(min_partitions)
        vds = hl.vds.read_vds(vds_url, intervals=intervals)
        print("Placed command for re-partitioning and Hail will execute the repartitioning later as Hail is lazy.")
    else:
        print("No repartitioning required.")
    return vds


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


def perform_chr_filtering_if_needed(vds: hl.vds.VariantDataset, is_test: bool, test_chrs: str) -> hl.vds.VariantDataset:
    """
    Handle the filtering of the Variant Dataset (VDS) for optimization and test purposes.

    Parameters:
    vds (hl.vds.VariantDataset:): The input VDS.
    is_test (bool): Flag indicating if this is a test run. Ignored if is_test is False. A test run will only work on
                    chromosomes specified in the argument "test_chrs" in the workflow. It will work
                    on all chr1-Y chromosomes if it is not a test run.
    test_chrs (Optional[str]): Specific chromosome to filter, if applicable.

    Returns:
    hl.vds.VariantDataset: The processed and optionally repartitioned VDS.
    """
    if is_test:
        test_chrs_list = re.split(r'[,_\s]+', test_chrs)
        chrs_grch38 = hl.get_reference('GRCh38').contigs
        chrs_grch38 = chrs_grch38[0:24]
        validate_input_chromosomes(test_chrs_list, chrs_grch38)
        vds = hl.vds.filter_chromosomes(vds, keep=test_chrs_list)
        print(f"Filtering for chromosome: {test_chrs}")
    else:
        print("No chromosome filtering applied.")
    return vds


def annotate_beds_and_filter(vds: hl.vds.VariantDataset, bed_file_dict: Dict[str, str]) -> hl.vds.VariantDataset:
    """
    Annotate and filter the Variant Dataset (VDS) based on regions specified in UCSC BED files.

    Parameters:
    vds (hl.vds.VariantDataset): The input Variant Dataset.
    bed_file_dict (Dict[str, str]): A dictionary mapping types of UCSC BED files to their file paths.

    Returns:
    hl.vds.VariantDataset: The annotated and filtered Variant Dataset.
    """

    def join_bed(vds: hl.vds.VariantDataset, bed_path: str, annotation_name: str) -> hl.vds.VariantDataset:
        """
        Parameters:
        vds (hl.vds.VariantDataset): The Hail Variant Dataset.
        bed_path (str): The gs URL to the UCSC bed file.
        annotation_name: Annotation name given to the UCSC bed file while annotating the intervals in
                         the bed file to the VDS.
        Returns:
        hl.vds.VariantDataset: An annotated Variant Dataset.
        """
        ht = hl.import_bed(bed_path, reference_genome=vds.reference_data.locus.dtype.reference_genome)
        print(f'UCSC Bed file at {bed_path} contains {ht.count()} intervals')
        vd = vds.variant_data
        vd = vd.annotate_rows(**{annotation_name: hl.is_defined(ht[vd.locus])})
        return hl.vds.VariantDataset(reference_data=vds.reference_data, variant_data=vd)

    for bed_name, bed_path in bed_file_dict.items():
        vds = join_bed(vds, bed_path, bed_name)

    beds = bed_file_dict.keys()
    vds = hl.vds.VariantDataset(
        reference_data=vds.reference_data,
        variant_data=vds.variant_data.filter_rows(hl.any(*(vds.variant_data[fd] for fd in beds)), keep=True)
    )
    return vds


# this will be extracted into a shared library in the future
def make_dense_mt(vds: hl.vds.VariantDataset, is_filtering_FT: bool, is_keep_as_vqsr_fields: bool,
                  max_alt_alleles: Optional[int], fields_to_drop: str) -> hl.MatrixTable:
    """
    Convert a Variant Dataset (VDS) to a densified MatrixTable (MT) with additional annotations, and drop some
    fields if fields exist in the VDS.

    Parameters:
    vds (hl.vds.VariantDataset): The input Variant Dataset.
    is_filtering_FT (bool): Flag to filter GTs based on the FT field. If "FT" is "Fail", then set
                        corresponding "GT" to missing.
    is_keep_as_vqsr_fields (bool): Flag to keep VQSR fields.
    max_alt_alleles (Optional[int]): Maximum number of alternate alleles for filtering.
    fields_to_drop (str): Specific fields to drop from the VDS while converting to a dense MT, e.g., "as_vqsr, LAD, tranche_data,
                          truth_sensitivity_snp_threshold, truth_sensitivity_indel_threshold,
                          snp_vqslod_threshold, indel_vqslod_threshold", delimiter can be , _ or space.

    Returns:
    hl.MatrixTable: The resulting densified and annotated MatrixTable.
    """
    vd_gt = vds.variant_data

    if max_alt_alleles is not None:
        vd_gt = vd_gt.filter_rows(hl.len(vd_gt.alleles) < max_alt_alleles)

    vd_gt = vd_gt.annotate_entries(AD=hl.vds.local_to_global(vd_gt.LAD, vd_gt.LA, n_alleles=hl.len(vd_gt.alleles), fill_value=0, number='R'))
    vd_gt = vd_gt.transmute_entries(GT=hl.vds.lgt_to_gt(vd_gt.LGT, vd_gt.LA))

    if 'FT' in vd_gt.entry:
        vd_gt = vd_gt.transmute_entries(FT=hl.if_else(vd_gt.FT, "PASS", "FAIL"))

    if 'gvcf_info' in vd_gt.entry:
        vd_gt = vd_gt.drop('gvcf_info')

    d_callset = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, vd_gt))
    
    if is_filtering_FT and "FT" in d_callset.entry:
        d_callset = d_callset.annotate_entries(GT=hl.or_missing((~hl.is_defined(d_callset.FT)) | (d_callset.FT.contains("PASS")), d_callset.GT))

    d_callset = hl.variant_qc(d_callset)
    d_callset = d_callset.annotate_rows(info=hl.struct(AC=d_callset.variant_qc.AC[1:], AF=d_callset.variant_qc.AF[1:], AN=d_callset.variant_qc.AN, homozygote_count=d_callset.variant_qc.homozygote_count))

    if is_keep_as_vqsr_fields and ('as_vqsr' in d_callset.row):
        d_callset = d_callset.annotate_rows(info=d_callset.info.annotate(AS_VQSLOD=d_callset.as_vqsr.values().vqslod, AS_YNG=d_callset.as_vqsr.values().yng_status))

    d_callset = d_callset.annotate_rows(variant_qc=d_callset.variant_qc.drop("AC", "AF", "AN", "homozygote_count"))
    fields_to_drop_list = re.split(r'[,_\s]+', fields_to_drop)
    d_callset = d_callset.drop(*(f for f in fields_to_drop_list if f in d_callset.entry or f in d_callset.row or f in d_callset.col))

    return d_callset


# Command-line argument parsing
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process VDS and UCSC BED files to create a filtered MatrixTable.")
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
    parser.add_argument('--test_chrs', help='Chromosome to filter, e.g., "chr22_chr21"', required=False)
    parser.add_argument('--version_label', help='Version identifier for the output', required=True)
    parser.add_argument('--task_label', help='task performed in the script, e.g.: "vds2mt"', required=True)
    parser.add_argument('--min_partitions', type=int, help='Minimum number of partitions for parallelization', required=True)
    parser.add_argument('--vds_url', help='URL of the input VDS', required=True)
    parser.add_argument('--is_filtering_FT', type=lambda x: (str(x).lower() == 'true'), help='Filter GTs based on "Failed" FT field', required=True)
    parser.add_argument('--is_keep_as_vqsr_fields', type=lambda x: (str(x).lower() == 'true'), help='Keep VQSR fields', required=True)
    parser.add_argument('--max_alt_alleles', type=int, help='Maximum number of alternate alleles for filtering', required=False)
    parser.add_argument('--staging_bucket', help='Temporary bucket to hold files on dataproc cluster', required=True)
    parser.add_argument('--basis_mt_json_file_name', help='File name of the basis dense MatrixTables json file, '
                                                          'suffix should be .json', required=True)
    parser.add_argument('--fields_to_drop', help='Fields to drop from the VDS while converting to a '
                                                 'dense MT, e.g., "as_vqsr, LAD, tranche_data, '
                                                 'truth_sensitivity_snp_threshold, '
                                                 'truth_sensitivity_indel_threshold, '
                                                 'snp_vqslod_threshold, indel_vqslod_threshold", '
                                                 'delimiter can be , _ or space.', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Initialize Hail with Spark settings
    # The Hail team will likely make this the default behavior in a later version.
    hail_init(args.executor_memory, args.executor_cores, args.driver_memory, args.reference_genome)

    # Process VDS for repartitioning if needed
    vds_processed = perform_repartition_if_needed(args.min_partitions, args.vds_url)

    # Process VDS for chromosome filtering if it's a test run
    vds_processed = perform_chr_filtering_if_needed(vds_processed, args.is_test, args.test_chrs)

    # Create a dictionary for UCSC BED files { "exome": exome_bed_path, "clinvar": clinvar_bed_path,.. etc. }
    bed_file_dict = create_bed_dict(args.exome_bed_path, args.clinvar_bed_path, args.acaf_threshold_bed_path)

    # Annotate the processed VDS with UCSC BED file regions and filter it
    vds_annotated = annotate_beds_and_filter(vds_processed, bed_file_dict)

    # Convert the annotated VDS to a densified MatrixTable with additional annotations
    mt = make_dense_mt(vds_annotated, args.is_filtering_FT, args.is_keep_as_vqsr_fields, args.max_alt_alleles, args.fields_to_drop)

    # Generate the output path for the Hail MatrixTable
    basis_mt_dict = create_output_path_dict(args.is_test, args.test_chrs, args.task_label, args.version_label, args.output_gs_url)

    # Write basis mt dictionary to a json file as output
    write_dict_to_json(basis_mt_dict, args.basis_mt_json_file_name, args.output_gs_url,
                       args.task_label, args.version_label)

    # Export the resulting MatrixTable to the specified location
    mt.write(basis_mt_dict["dense_mt"], overwrite = True)

