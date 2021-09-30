#!/usr/bin/env python3

import argparse
import gzip
import numpy as np

import scipy.io
import scipy.sparse


def main():
    description = """Create npz, npy file from the mtx files produced by STARsolo"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--barcodes',
                        dest='barcodes',
                        required=True,
                        help="The barcodes file")

    parser.add_argument('--features',
                        dest='features',
                        required=True,
                        help="The features file")

    parser.add_argument('--matrix',
                        dest='matrix',
                        required=True,
                        help="The matrix file")

    args = parser.parse_args()

    # covert the mtx file to the matrix
    matrix = scipy.io.mmread(args.matrix).transpose().tocsr()
    nonzero_row_indices, _ = matrix.nonzero()
    unique_nonzero_row_indices = np.sort(np.unique(nonzero_row_indices))
    # we need to keep only those rows that have non-zero reads/counts
    scipy.sparse.save_npz("sparse_counts.npz", matrix[unique_nonzero_row_indices, :], compressed=True)

    # read the barcodes file and create the barcode to index
    barcodes = []
    with gzip.open(args.barcodes, 'rt') if args.barcodes.endswith('.gz') else  \
        open(args.barcodes, 'r') as fin:
        for line in fin:
           if line.startswith(r'^#'): # skip comments
              continue
           fields = line.strip().split('\t')
           barcodes.append(fields[0])

    row_index = np.asarray(barcodes)
    # we need to keep only those barcodes that have non-zero reads/counts
    np.save("sparse_counts_row_index.npy", row_index[unique_nonzero_row_indices])
      
    # read the features file and create the feature to index map
    features = []
    with gzip.open(args.features, 'rt') if args.features.endswith('.gz') else  \
        open(args.features, 'r') as fin:
        for line in fin:
           if line.startswith(r'^#'): # skip comments
              continue
           fields = line.strip().split('\t')
           features.append(fields[0])

    row_index = np.asarray(features)
    np.save("sparse_counts_col_index.npy", row_index)

if __name__ == '__main__':
    main()
