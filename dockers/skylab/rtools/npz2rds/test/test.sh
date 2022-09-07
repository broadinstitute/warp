#!/bin/bash

set -eo pipefail

## The location of the executatbles
execDir=..

## Add the exec directory to the PATH
PATH=$execDir:$PATH

## Unzip the input test data
printf "Extracting test data... "
tar xzf input.tar.gz
printf "done\n"

printf "Running Command...\n"
## Run the command
$execDir/npz2rds.sh -c input/sparse_counts_col_index.npy \
		    -r input/sparse_counts_row_index.npy \
		    -d input/sparse_counts.npz \
		    -o counts.rds
printf "done\n"

printf "Validating output checksum..."
## Check the output checksum
chksum=`md5sum counts.rds | cut -f 1 -d ' ' `
exitCode=0
if [[ "$chksum" == "f19b93dc8ecb651dd8956c44a2a53725" ]]
then
    printf "PASSED\n"
    exitCode=0
else
    printf "FAILED\n"
    exitCode=1
fi

## Cleanup
printf "Cleaning up..."
rm -r input/
rm counts.rds
printf "done\n"

exit $exitCode
