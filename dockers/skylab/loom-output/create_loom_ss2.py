import argparse
import csv
import os
import numpy as np
import scipy as sc
import loompy

def generate_col_attr(args):
    """Converts the QC of Smart Seq2 gene file pipeline outputs to loom file
    Args:
        qc_path (str): path to the QCs csv
    """
    # read the QC values
    qc_path = [p for p in args.qc_files if p.endswith("_QCs.csv")][0]
    with open(qc_path, 'r') as f:
        qc_values = [row for row in csv.reader(f)]

    metadata_labels = qc_values[0][1:]
    metadata_values = qc_values[2][1:]
    cell_id = qc_values[2][0]

    string_metadata = {}
    numeric_metadata = {}

    for label, value in zip(metadata_labels, metadata_values):

        # See if this is numeric or string
        numeric_value = None
        try:
            numeric_value = float(value)
        except ValueError:
            try:
                numeric_value = float(value.strip("%"))/100
            except ValueError:
                pass

        if numeric_value is not None:
            numeric_metadata[label] = numeric_value
        else:
            string_metadata[label] = value

    # Metrics
    # Write the string and numeric metadata separately
    sorted_string_labels = sorted(string_metadata.keys())
    sorted_string_values = [string_metadata[m] for m in sorted_string_labels]
    sorted_numeric_labels = sorted(numeric_metadata.keys())
    sorted_numeric_values = [numeric_metadata[m] for m in sorted_numeric_labels]

    # Column attributes
    col_attrs = dict()
    col_attrs["cell_names"] = [cell_id]
    col_attrs["CellID"] = [cell_id]
    col_attrs['input_id'] = [args.input_id]
    if args.input_name is not None:
        col_attrs['input_name'] = [args.input_name]

    if args.input_id_metadata_field:
        col_attrs["input_id_metadata_field"] = [args.input_id_metadata_field]
    if args.input_name_metadata_field:
            col_attrs["input_name_metadata_field"] = [args.input_name_metadata_field]

    numeric_field_names = np.array(sorted_numeric_labels[:])
    for i in range(0, numeric_field_names.shape[0]):
        name = numeric_field_names[i]
        data = np.array([sorted_numeric_values])[:,i]
        col_attrs[name] = data
    string_field_names = np.array(sorted_string_labels)
    for i in range(0, string_field_names.shape[0]):
        name = string_field_names[i]
        data = np.array([sorted_string_values])[:,i]
        col_attrs[name] = data


    return col_attrs


def generate_csr_spase_coo(expr):

    nrows, ncols = np.shape(expr)
    expr_coo = sc.sparse.coo_matrix(expr[:])
    xcoord = []
    ycoord = []
    value = []

    for k in range(0, expr_coo.data.shape[0]):
        xcoord.append(expr_coo.row[k])
        ycoord.append(expr_coo.col[k])
        value.append(expr_coo.data[k])

    xcoord = np.asarray(xcoord)
    ycoord = np.asarray(ycoord)
    value = np.asarray(value)

    expr_sp_t = sc.sparse.coo_matrix((value, (ycoord, xcoord)), shape=(ncols, nrows))
    return expr_sp_t


def generate_row_attr_and_matrix(rsem_gene_results_path):
    """Converts the Smart Seq2 gene file pipeline outputs to loom file
    Args:
        input_path (str): file where the SS2 pipeline expression counts are
    """
    reader = csv.DictReader(open(rsem_gene_results_path), delimiter="\t")
    expression_values = {}
    count_values = {}
    for row in reader:
        expression_values[row["gene_id"]] = float(row["TPM"])
        count_values[row["gene_id"]] = float(row["expected_count"])
    sorted_gene_ids = sorted(expression_values.keys())
    sorted_tpms = [expression_values[g] for g in sorted_gene_ids]
    sorted_counts = [count_values[g] for g in sorted_gene_ids]
    expression_tpms = generate_csr_spase_coo( [sorted_tpms])
    expected_counts = generate_csr_spase_coo([sorted_counts])
    row_attrs = {"ensembl_ids":np.array(sorted_gene_ids),"gene_names":np.array(sorted_gene_ids),"Gene":np.array(sorted_gene_ids)}

    return row_attrs, expression_tpms,expected_counts


def create_loom_files(args):
    """This function creates the loom file or folder structure in output_loom_path in
       format file_format, with input_id from the input folder analysis_output_path
    Args:
        input_id (str): sample or cell id
        qc_analysis_output_files_string (str): a string with the file names in the QCGroup of SS2
            pipeline output, separated by commas
        rsem_genes_results_file (str): the file for the expression count
        output_loom_path (str): location of the output loom
    """
    # generate a dictionary of column attributes
    col_attrs =  generate_col_attr(args)

    # add the expression count matrix data
    # generate a dictionary of row attributes
    row_attrs, expr_tpms, expr_counts = generate_row_attr_and_matrix(args.rsem_genes_results_file)

    attrDict = dict()
    attrDict['pipeline_version'] = args.pipeline_version

    #generate loom file
    loompy.create(args.output_loom_path, expr_tpms, row_attrs, col_attrs, file_attrs=attrDict)
    ds = loompy.connect(args.output_loom_path)
    ds.layers['estimated_counts'] = expr_counts
    ds.close()


def main():
    description = """This script converts the some of the SmartSeq2 pipeline outputs in to
                   loom format (https://linnarssonlab.org/loompy/format/index.html) relevant output.
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--qc_files',
                        dest="qc_files",
                        nargs = "+",
                        help=('the grouped QC files from the GroupQCOutputs task of SS2 '
                              'Single Sample workflow'))

    parser.add_argument('--rsem_genes_results',
                        dest="rsem_genes_results_file",
                        help='path to the folder containing the files to be added to the loom')

    parser.add_argument('--output_loom_path',
                        dest="output_loom_path",
                        help='path where the loom file is to be created')

    parser.add_argument('--input_id',
                        dest="input_id",
                        help='the sample name in the bundle')

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

    parser.add_argument('--pipeline_version',
                        default="Unknown sample",
                        help='the version of SS2 used to generate data')

    args = parser.parse_args()

    create_loom_files(args)


if __name__ == '__main__':
    main()
