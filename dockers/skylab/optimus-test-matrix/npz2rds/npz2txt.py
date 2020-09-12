#!/usr/bin/env python3

## Author: Nick Barkas <nbarkas@broadinstitute.org>
## Date: December 13, 2018
## Description: A python script that will read the output count matrix from
##   the Optimus pipeline and save it as discrete text files in the selected output
##   directory

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Convert Optimus count matrix to rds file")
parser.add_argument('--col-index',dest='colindex',help='column index npy file')
parser.add_argument('--row-index',dest='rowindex', help='row index npy file')
parser.add_argument('--counts',dest='counts', help='counts npz file')
parser.add_argument('--output-dir',dest='output',help='output directory')

args = parser.parse_args()

## TODO: Check that the args are valid

## Column indices
col_indices = np.load(args.colindex)
np.savetxt(args.output + '/sparse_counts_col_index.txt', col_indices, fmt='%s')

## Row indices
row_indices = np.load(args.rowindex)
np.savetxt(args.output + '/sparse_counts_row_index.txt', row_indices, fmt='%s')

## The count data
data = np.load(args.counts)

## Save each component as individual txt file
np.savetxt(args.output + '/sparse_counts_indices.txt',data['indices'],fmt='%i')
np.savetxt(args.output + '/sparse_counts_indptr.txt',data['indptr'],fmt='%i')
np.savetxt(args.output + '/sparse_counts_shape.txt',data['shape'],fmt='%i')
np.savetxt(args.output + '/sparse_counts_data.txt',data['data'],fmt='%i')

### End of Script
