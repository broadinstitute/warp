import argparse
import csv
import gzip
import re
import numpy as np
import loompy
from scipy import sparse
import pandas as pd
import scipy as sc
import logging


def create_gene_id_name_map(gtf_file):
    """ Creates a map from gene_id to gene_name by reading in the GTF file

    Args:
        gtf_file (str): annotation file

    Return:
        gene_id_name_map (Dict[str, str]): dictonary gene ids to gene names
    """
    gene_id_name_map = {}

    # loop through the lines and find the gene_id and gene_name pairs
    with gzip.open(gtf_file, "rt") if gtf_file.endswith(".gz") else open(
            gtf_file, "r"
    ) as fpin:
        for _line in fpin:
            line = _line.strip()
            gene_id_res = re.search(r"gene_id ([^;]*);", line)
            gene_name_res = re.search(r"gene_name ([^;]*);", line)

            if gene_id_res and gene_name_res:
                gene_id = gene_id_res.group(1).replace('"', "")
                gene_name = gene_name_res.group(1).replace('"', "")
                gene_id_name_map[gene_id] = gene_name

    return gene_id_name_map

def  generate_row_attr(args):
    """Converts the gene metrics from the Optimus pipeline to row attributes for loom file

    Args:
        verbose (bool): whether to output verbose messages for debugging purposes
    """

    # read the gene metrics names and values
    input_path = args.gene_metrics
    if input_path.endswith(".gz"):
        with gzip.open(input_path, "rt") as f:
            gene_metrics = [row for row in csv.reader(f)]
    else:
        with open(input_path, "r") as f:
            gene_metrics = [row for row in csv.reader(f)]

    if args.verbose:
        logging.info("# gene numeric metadata {0}".format(len(gene_metrics[0][1:])))

    # read the gene ids  and adds into the gene_ids dataset
    gene_ids = np.load(args.gene_id_1)
    # add the gene names for the gene ids
    if args.annotation_file and len(gene_ids):
        gene_id_name_map = create_gene_id_name_map(args.annotation_file)
        gene_names = [
            gene_id_name_map.get(gene_id, "") for gene_id in gene_ids
        ]
    else: print("error in gene_id and annotation file!")


    # Gene metric values, the row and column sizes
    gene_ids_location = {
        gene_id: index for index, gene_id in enumerate(gene_ids)
    }

    # ignore the first line with the metric names in text
    ncols = 0
    gene_id_to_metric_values = {}
    for row in gene_metrics:
        # only consider genes that are in the count matrix
        if not row[0] in gene_ids_location:
            continue

        row_values = []
        for value_string in row[1:]:
            # some of the standard deviation values do not exist for one reads matches
            try:
                value = np.float32(value_string)
            except ValueError:
                value = np.nan
            row_values.append(value)
        gene_id_to_metric_values[row[0]] = row_values

        # note that all of these lengths are assumed to be equal and this check is already done in the pipeline
        ncols = len(row_values)

    # now insert the metrics of the cells that are in count matrix, i.e.,  the global variable "cell_ids"

    gene_metric_values = []
    for gene_id in gene_ids:
        if gene_id in gene_id_to_metric_values:
            gene_metric_values.append(gene_id_to_metric_values[gene_id])
        else:
            # if no metrics for a cell present in the count matrix then fill them with np.nans
            gene_metric_values.append([np.nan] * ncols)

    nrows = len(gene_ids)

    # Generate row attribute dictionary
    row_attrs = { "gene_names":gene_names, "ensembl_ids": gene_ids, "Gene": gene_names}

    gene_metrics_data =np.array(gene_metric_values)
    numeric_field_names = gene_metrics[0][1:]
    for i in range(0, len(numeric_field_names)):
        name = numeric_field_names[i]
        data = gene_metrics_data[:, i]
        row_attrs[name] = data
    if args.verbose:
        logging.info("# of genes: {0}".format(nrows))
        logging.info("# of gene metadate metrics: {0}".format(ncols))
    return row_attrs

def generate_col_attr(barcode_1_path,barcode_2_path,cell_metrics_path):
    """Converts cell metrics from the Optimus pipeline to loom file

    Args:
        input_path (str): file containing gene metrics name and values
        cell_ids (list): list of cell ids
        verbose (bool): whether to output verbose messages for debugging purposes
    """
    barcode_1 = np.load(barcode_1_path)
    barcode_2 = np.load(barcode_2_path)
    cell_ids = np.union1d(barcode_1,barcode_2)

    # Read the csv input files
    metrics_df = pd.read_csv(cell_metrics_path, header=0,  dtype=str)
    # Check that input is valid
    if metrics_df.shape[0] == 0 or metrics_df.shape[1] == 0:
        logging.error("Cell metrics table is not valid")
        raise ValueError()
    metrics_df = metrics_df.rename(columns={"Unnamed: 0": "cell_id"})

    # Order the cells by merging with cell_ids
    cellorder_df = pd.DataFrame(data={"cell_id": cell_ids})

    # Split the pandas DataFrame into different data types for storing in the ZARR
    IntColumnNames = [  # UInt
        "n_reads",
        "noise_reads",
        "perfect_molecule_barcodes",
        "reads_mapped_exonic",
        "reads_mapped_intronic",
        "reads_mapped_utr",
        "reads_mapped_uniquely",
        "reads_mapped_multiple",
        "duplicate_reads",
        "spliced_reads",
        "antisense_reads",
        "n_molecules",
        "n_fragments",
        "fragments_with_single_read_evidence",
        "molecules_with_single_read_evidence",
        "perfect_cell_barcodes",
        "reads_mapped_intergenic",
        "reads_unmapped",
        "reads_mapped_too_many_loci",
        "n_genes",
        "genes_detected_multiple_observations"
    ]

    FloatColumnNames = [ # Float32
        "molecule_barcode_fraction_bases_above_30_mean",
        "molecule_barcode_fraction_bases_above_30_variance",
        "genomic_reads_fraction_bases_quality_above_30_mean",
        "genomic_reads_fraction_bases_quality_above_30_variance",
        "genomic_read_quality_mean",
        "genomic_read_quality_variance",
        "reads_per_fragment",
        "fragments_per_molecule",
        "cell_barcode_fraction_bases_above_30_mean",
        "cell_barcode_fraction_bases_above_30_variance",
        "n_mitochondrial_genes",
        "n_mitochondrial_molecules",
        "pct_mitochondrial_molecules"
    ]

    newcolnames = []

    final_df = cellorder_df.merge(metrics_df, on="cell_id", how="left")
    ColumnNames = IntColumnNames + FloatColumnNames
    BoolColumnNames = []

    # Split the dataframe
    final_df_non_boolean = final_df[ColumnNames]
    # Create a numpy array for the column names
    final_df_bool_column_names = final_df[BoolColumnNames].columns.values
    # Create a numpy array of the same shape of booleans
    final_df_bool = np.full(final_df[BoolColumnNames].shape, True)


    final_df_non_boolean = final_df_non_boolean.apply(pd.to_numeric)

    # COLUMN/CELL Metadata
    col_attrs = dict()
    col_attrs["cell_names"] = cell_ids
    col_attrs["CellID"] = cell_ids
    bool_field_names = final_df_bool_column_names

    # Create metadata tables and their headers for bool
    for i in range(0, bool_field_names.shape[0]):
        name = bool_field_names[i]
        data = final_df_bool[:, i]
        col_attrs[name] = data

    # Create metadata tables and their headers for float
    float_field_names = list(final_df_non_boolean.columns)

    for i in range(len(float_field_names)):
        name = float_field_names[i]
        data = final_df_non_boolean[name].to_numpy()
        col_attrs[name] = data

    return col_attrs

def generate_matrix(count_matrix_path):

    # read .npz file expression counts and add it to the expression_counts dataset
    exp_counts = np.load(count_matrix_path)
    # now convert it back to a csr_matrix object
    csr_exp_counts = sparse.csr_matrix(
        (exp_counts["data"], exp_counts["indices"], exp_counts["indptr"]),
        shape=exp_counts["shape"],
    )

    nrows, ncols = csr_exp_counts.shape
    expr_sp = sc.sparse.coo_matrix((nrows, ncols), np.float32)

    xcoord = []
    ycoord = []
    value = []

    chunk_row_size = 10000
    chunk_col_size = 10000

    for i in range(0, nrows, chunk_row_size):
        for j in range(0, ncols, chunk_col_size):
            p = chunk_row_size
            if i + chunk_row_size > nrows:
                p = nrows - 1
            q = chunk_col_size
            if j + chunk_col_size > ncols:
                q = ncols - j
            expr_sub_row_coo = sc.sparse.coo_matrix(csr_exp_counts[i:i + p, j:j + q].toarray())
            for k in range(0, expr_sub_row_coo.data.shape[0]):
                xcoord.append(expr_sub_row_coo.row[k] + i)
                ycoord.append(expr_sub_row_coo.col[k] + j)
                value.append(expr_sub_row_coo.data[k])

    xcoord = np.asarray(xcoord)
    ycoord = np.asarray(ycoord)
    value = np.asarray(value)

    expr_sp_t = sc.sparse.coo_matrix((value, (ycoord, xcoord)), shape=(expr_sp.shape[1], expr_sp.shape[0]))

    del xcoord
    del ycoord
    del value

    return expr_sp_t

# This function takes two sparse matrices with different number of columns and makes them the same shape
def changeCoord(sp_1,sp_2, barcodes_1, barcodes_2):
    barcodes = np.union1d(barcodes_1,barcodes_2)
    i1=np.where(np.in1d(barcodes,barcodes_1))[0]
    i2= np.where(np.in1d(barcodes,barcodes_2))[0]
    sp_1_col_new = [i1[i] for i in sp_1.col]
    sp_2_col_new = [i2[i] for i in sp_2.col]
    sp_1_new = sc.sparse.coo_matrix((sp_1.data, (sp_1.row, sp_1_col_new)), shape=(sp_1.shape[0],len(barcodes)))
    sp_2_new = sc.sparse.coo_matrix((sp_2.data, (sp_2.row, sp_2_col_new)), shape=(sp_2.shape[0],len(barcodes)))

    return sp_1_new, sp_2_new

def create_loom_files(args):
    """This function creates the loom file or folder structure in output_loom_path in format file_format,
       with input_id from the input folder analysis_output_path

    Args:
        args (argparse.Namespace): input arguments for the run
    """
    version = "1.0.0"


    # generate a dictionary of row attributes
    row_attrs =  generate_row_attr(args)

    # generate a dictionarty of column attributes
    col_attrs =  generate_col_attr(args.cell_id_1,args.cell_id_2,args.cell_metrics)

    # add the expression count matrix data
    sp1 = generate_matrix(args.count_matrix_1)
    sp2 = generate_matrix(args.count_matrix_2)

    barcode_1 = np.load(args.cell_id_1)
    barcode_2 = np.load(args.cell_id_2)
    expr_sp_t, exon_sp_t = changeCoord(sp1,sp2,barcode_1,barcode_2)

# add input_id to col_attrs
    col_attrs['input_id'] = np.repeat(args.input_id, expr_sp_t.shape[1])

    # generate global attributes
    attrDict = dict()
    attrDict['expression_data_type'] = args.expression_data_type
    attrDict['optimus_output_schema_version'] = version
    attrDict['input_id'] = args.input_id
    if args.input_name is not None:
        attrDict['input_name'] = args.input_name
    if args.input_id_metadata_field is not None:
        attrDict['input_id_metadata_field'] = args.input_id_metadata_field
    if args.input_name_metadata_field is not None:
        attrDict['input_name_metadata_field'] = args.input_name_metadata_field
    attrDict['pipeline_version'] = args.pipeline_version
    #generate loom file
    loompy.create(args.output_loom_path, expr_sp_t, row_attrs, col_attrs, file_attrs=attrDict)

    ds = loompy.connect(args.output_loom_path)
    ds["exon_counts"] = exon_sp_t

    ds.close()

def main():
    description = """This script converts the some of the Optimus outputs in to
                   loom format (http://linnarssonlab.org/loompy/index.html) relevant output.
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--cell_metrics",
        dest="cell_metrics",
        required=True,
        help="a .csv file path for the cell metrics, an output of the MergeCellMetrics task",
    )

    parser.add_argument(
        "--gene_metrics",
        dest="gene_metrics",
        required=True,
        help="a .csv file path for the gene metrics, an output of the MergeGeneMetrics task",
    )

    parser.add_argument(
        "--cell_id_1",
        dest="cell_id_1",
        required=True,
        help="a .npy file path for the cell barcodes in the main layer, an output of the MergeCountFiles task",
    )

    parser.add_argument(
        "--gene_id_1",
        dest="gene_id_1",
        required=True,
        help="a .npy file path for the gene ids in the main layer of the loom, an output of the MergeCountFiles task",
    )

    parser.add_argument(
        "--count_matrix_1",
        dest="count_matrix_1",
        required=True,
        help="a .npz file path for the count matrix for the main layer of the loom, an output of the MergeCountFiles task",
    )
    parser.add_argument(
        "--cell_id_2",
        dest="cell_id_2",
        required=True,
        help="a .npy file path for the cell barcodes for the second layer, an output of the MergeCountFiles task",
    )

    parser.add_argument(
        "--gene_id_2",
        dest="gene_id_2",
        required=True,
        help="a .npy file path for the gene ids for the second layer, an output of the MergeCountFiles task",
    )

    parser.add_argument(
        "--count_matrix_2",
        dest="count_matrix_2",
        required=True,
        help="a .npz file path for the count matrix for the second layer, an output of the MergeCountFiles task",
    )

    parser.add_argument(
        "--annotation_file",
        dest="annotation_file",
        default=None,
        required=False,
        help="annotation file in GTF format",
    )

    parser.add_argument(
        "--output_path_for_loom",
        dest="output_loom_path",
        required=True,
        help="path to .loom file is to be created",
    )

    parser.add_argument(
        "--input_id",
        dest="input_id",
        required=True,
        help="the sample name in the bundle",
    )

    parser.add_argument(
        "--input_name",
        dest="input_name",
        help= "sequencing_input.biomaterial_core.biomaterial_id in HCA metadata, defined by the user",
    )

    parser.add_argument(
        "--input_id_metadata_field",
        dest="input_id_metadata_field",
        help= "sequencing_process.provenance.document_id: [UUID] defined by the user",
    )

    parser.add_argument(
        "--input_name_metadata_field",
        dest="input_name_metadata_field",
        help= "sequencing_input.biomaterial_core.biomaterial_id defined by the user",
    )

    parser.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true",
        help="whether to output verbose debugging messages",
    )

    parser.add_argument(
        "--expression_data_type",
        dest="expression_data_type",
        default="exonic",
        choices=["exonic", "whole_transcript"],
        help="The expression data type",
    )

    parser.add_argument(
        "--pipeline_version",
        dest="pipeline_version",
        required=True,
        help="The version of Optimus that generated data",
    )

    args = parser.parse_args()

    create_loom_files(args)

if __name__ == "__main__":
    main()
