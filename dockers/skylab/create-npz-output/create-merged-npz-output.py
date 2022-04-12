#!/usr/bin/env python3

import argparse
import gzip
import numpy as np
import sys

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

def main():
    description = """Create npz, npy file from the mtx files produced by STARsolo"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--barcodes',
                        dest='barcodes',
                        nargs = '+',
                        required=True,
                        help="The barcodes file")

    parser.add_argument('--features',
                        dest='features',
                        nargs = '+',
                        required=True,
                        help="The features file")

    parser.add_argument('--matrix',
                        dest='matrix',
                        nargs = '+',
                        required=True,
                        help="The matrix file")

    args = parser.parse_args()
    if(len(args.barcodes) != len(args.matrix)):
      print("Number of barcode files are not the same as number of count matrices.")
      exit(1)
    elif(len(args.features) != len(args.matrix)):
      print("Number of feature files are not the same as number of count matrices.")
      exit(1)
    elif(len(args.barcodes) != len(args.features)):
      print("Number of barcode files are not the same as number of feature files.")
      exit(1)

    # features list file 
    features_list = []
    for i in range(len(args.features)):
      features_list = np.union1d(features_list, read_features_list(args.features[i]))
      print(len(features_list))
    np.save("sparse_counts_col_index.npy", features_list)
    # read the barcodes file and create the barcode to index
    barcodes_list = []
    for i in range(len(args.barcodes)):
      barcodes_list = np.union1d(barcodes_list, read_barcodes_list(args.barcodes[i]))
      print(len(barcodes_list))

    np.save("sparse_counts_row_index.npy", barcodes_list)
    data_dict = {}
    for i in range(len(args.matrix)):
        # convert the mtx file to the matrix
        matrix = scipy.io.mmread(args.matrix[i]).transpose().tocoo()

        for col, row, count in zip(matrix.col, matrix.row, matrix.data):
             if (row, col) not in data_dict:
                 data_dict[(row,col)] = 0
             data_dict[(row,col)] += count
        
    #unique_nonzero_row_indices = np.sort(np.unique(nonzero_row_indices))
    scipy.sparse.save_npz("sparse_counts.npz", matrix, compressed=True)


if __name__ == '__main__':
    main()
