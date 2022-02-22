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



def generate_matrix(args):

    # read .npz file expression counts and add it to the expression_counts dataset
    exp_counts = np.load(args.count_matrix)
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

def create_loom_files(args):
    """This function creates the loom file or folder structure in output_loom_path in format file_format,
       with input_id from the input folder analysis_output_path

    Args:
        args (argparse.Namespace): input arguments for the run
    """
    version = "1.0.0"

    # add the expression count matrix data
    expr_sp_t = generate_matrix(args)


    # generate global attributes

    #generate loom file
    loompy.create(args.output_loom_path, expr_sp_t, row_attrs, col_attrs, file_attrs=attrDict)

def main():
    description = """This script adds exon counts to the Optimus outputs for single nucleus mode in to
                   loom format (http://linnarssonlab.org/loompy/index.html) relevant output.
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--gene_metrics",
        dest="gene_metrics",
        required=True,
        help="a .csv file path for the gene metrics, an output of the MergeGeneMetrics task",
    )

    parser.add_argument(
        "--cell_id",
        dest="cell_ids",
        required=True,
        help="a .npy file path for the cell barcodes, an output of the MergeCountFiles task",
    )

    parser.add_argument(
        "--gene_id",
        dest="gene_ids",
        required=True,
        help="a .npy file path for the gene ids, an output of the MergeCountFiles task",
    )

    parser.add_argument(
        "--count_matrix",
        dest="count_matrix",
        required=True,
        help="a .npz file path for the count matrix, an output of the MergeCountFiles task",
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




    args = parser.parse_args()

    create_loom_files(args)

if __name__ == "__main__":
    main()
