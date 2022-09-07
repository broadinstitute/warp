#!/bin/bash

## Description: Unit test for emptyDropsWrapper script

## Test parameters
testDataURL="http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
testDataFileName="pbmc4k_raw_gene_bc_matrices.tar.gz"
testDataInputMatrixPath="raw_gene_bc_matrices/GRCh38"
rdsFileName="pbmc4k.rds"
emptyDropsOutput="pbmc4k_emptyDrops.csv"
md5checksum="3de6a25eae8d5522c3a6989ab92479c1"
# Extra things extracted from the "testDataURL" you want to clean up
extraCleanup="raw_gene_bc_matrices"

## Download some sample data from 10X
printf "Downloading pbmc4k data..."
wget -q $testDataURL
printf "done\n"

## Decompress dataset
printf "Decompressing archive..."
tar xzf $testDataFileName
printf "done\n"

## Prepare RDS files for reading from the the 10X matrix
printf "Converting mtx to rds..."
./prepRDS.R --input ${testDataInputMatrixPath} --output ${rdsFileName}
printf "done\n"

## Run empty drops
printf "Running emptyDrops..."
../emptyDropsWrapper.R -i ${rdsFileName} -o ${emptyDropsOutput} --emptydrops-lower 1
printf "done\n"

## Check the output md5 checksum
## Note that empty drops is based on Monte Carlo simulation and therefore
## the output is not guaranteed to be deterministic, in our
## test it was and therefore the following checksum always
## validates
printf "Verifying checksum..."
md5out=`md5sum ${emptyDropsOutput} | cut -f 1 -d ' '`

exitCode=0
if [ "$md5out" = "$md5checksum" ];
then
    echo 'PASSED'
    exitCode=0
else
    echo "FAILED Expected $md5checksum got $md5out"
    exitCode=1
fi
printf "done\n"

## Cleanup
#rm ${emptyDropsOutput} ${testDataFileName} ${rdsFileName}
#rm -r ${extraCleanup}

exit $exitCode
