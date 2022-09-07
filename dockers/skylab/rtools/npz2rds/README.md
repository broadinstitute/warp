# npz2rds 
Converts Optimus output to an output suitable for reading in with R (rds file). Uses two scripts to do the conversion via txt files.

## Requirements

You need the following R packages installed: 'optparse'

The utility has been tested with the following:

R >= 3.5.1
python >= 3.6.2
- R >= 3.5.1
- python >= 3.6.2
- numpy >= 1.15.4

It probably works with much earlier versions.

On a mac you need the following for the test suite to work:

``
brew install md5sha1sum
``

# Example Usage
``
npz2rds.sh -c input/sparse_counts_col_index.npy -r input/sparse_counts_row_index.npy -d input/sparse_counts.npz -o counts.rds
``
