#!/usr/bin/env python3

import argparse
import gzip
import numpy as np
import sys
import scipy
import scipy.io
import scipy.sparse


def read_barcodes_list(barcodes_file):
    # read the barcodes file and create the barcode to index
    barcodes = []
    with gzip.open(barcodes_file, 'rt') if barcodes_file.endswith('.gz') else  \
        open(barcodes_file, 'r') as fin:
        for line in fin:
           if line.startswith(r'^#'): # skip comments
              continue
           fields = line.strip().split('\t')
           barcodes.append(fields[0])
    barcode_list = np.asarray(barcodes)
    return barcode_list

def read_features_list(features_file):
    # read the features file and create the feature to index map
    features = []
    with gzip.open(features_file, 'rt') if features_file.endswith('.gz') else  \
        open(features_file, 'r') as fin:
        for line in fin:
           if line.startswith(r'^#'): # skip comments
              continue
           fields = line.strip().split('\t')
           features.append(fields[0])
    features_list = np.asarray(features)
    return features_list

def changeCoord(feature_subset, features_union,
                barcode_subset, barcode_union,
                matrix):
    row_index = np.where(np.in1d(features_union,feature_subset))[0]
    col_index = np.where(np.in1d(barcode_union,barcode_subset))[0]
    mapped_row = []
    mapped_col = []
    data = []
    for r, c, v in zip(matrix.row, matrix.col,matrix.data):
        mapped_row.append(row_index[r])
        mapped_col.append(col_index[c])
        data.append(v)
    sp = scipy.sparse.coo_matrix((data, (mapped_row, mapped_col)), shape=(len(features_union),len(barcode_union)))
    return sp.transpose()

def main():
    description = """Create npz, npy file from the mtx files produced by STARsolo"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--barcodes',
                        dest='barcodes',
                        nargs='+',
                        required=True,
                        help="The barcodes file")

    parser.add_argument('--features',
                        dest='features',
                        nargs='+',
                        required=True,
                        help="The features file")

    parser.add_argument('--matrix',
                        dest='matrix',
                        nargs='+',
                        required=True,
                        help="The matrix file")
    parser.add_argument(
        "--input_id",
        dest="input_id",
        required=True,
        help="the sample name in the bundle",
    )

    args = parser.parse_args()
    if len(args.barcodes) != len(args.matrix):
        print("Number of barcode files are not the same as number of count matrices.")
        exit(1)
    elif len(args.features) != len(args.matrix):
        print("Number of feature files are not the same as number of count matrices.")
        exit(1)
    elif len(args.barcodes) != len(args.features):
        print("Number of barcode files are not the same as number of feature files.")
        exit(1)

    # features list file 
    features_dict = {}
    for i in range(len(args.features)):
        features_dict["sample{0}".format(i)] = read_features_list(args.features[i])
    features_list = read_features_list(args.features[0])

    # read the barcodes file and create the barcode to index
    barcodes_dict = {}
    for i in range(len(args.barcodes)):
        barcodes_dict["sample{0}".format(i)] = read_barcodes_list(args.barcodes[i])
    barcodes_list = []
    for key in barcodes_dict.keys():
        barcodes_list = np.union1d(barcodes_list, barcodes_dict[key])

    matrix_dict = {}
    for i in range(len(args.matrix)):
        matrix_dict["sample{0}".format(i)] = scipy.io.mmread(args.matrix[i])
    expr_sp = scipy.sparse.coo_matrix((len(barcodes_list),len(features_list)), np.float32)
    for k in features_dict.keys():
        sp = changeCoord(features_dict[k],features_list,
                         barcodes_dict[k],barcodes_list,
                         matrix_dict[k])
        expr_sp = expr_sp+sp

    matrix = expr_sp.tocsr()
    prev_index = 0
    discard_rows_indices = []
    for i in range(len(matrix.indptr)-1):
        if matrix.indptr[i] == matrix.indptr[i+1]:
            discard_rows_indices.append(i)

    nonzero_barcodes = [None] * (len(barcodes_list)-len(discard_rows_indices))
    j = 0
    i_dst = 0
    for i_src in range(len(barcodes_list)):
        if (j >= len(discard_rows_indices)) or (discard_rows_indices[j] != i_src):
            nonzero_barcodes[i_dst] = barcodes_list[i_src]
            i_dst = i_dst+1
        else:
            j = j+1

    reshaped_matrix = scipy.sparse.csr_matrix((matrix.data, matrix.indices, np.unique(matrix.indptr)),
                                              shape=(len(nonzero_barcodes), len(features_list)))

    scipy.sparse.save_npz(args.input_id+"_sparse_counts.npz", reshaped_matrix, compressed=True)
    np.save(args.input_id+"_sparse_counts_col_index.npy", features_list)
    np.save(args.input_id+"_sparse_counts_row_index.npy", nonzero_barcodes)

if __name__ == '__main__':
    main()
