#!/bin/bash

## Author: Nick Barkas <nbarkas@broadinstitute.org>
## Date: December 12, 2018
## Description: A wrapper script that converts optimus
##   output to RDS R files

set -eo pipefail

show_help() {
    echo "Usage: $0 [arguments]"
    echo ""
    echo "Arguments:"
    echo "  -c             column name file (.npy)"
    echo "  -r             row name file (.npy)"
    echo "  -d             data file (.npz)"
    echo "  -o             output file (.rds)"
    echo "  -t             optional temporary directory"
    echo "  -h             print this helpful message"
    echo ""
}


## Process command line arguments
OPTIND=1

## Init variables
colindexfile=""
rowindexfile=""
countsfile=""
outputfile=""
tmpdir=""

while getopts "hc:r:d:o:t:" opt; do
    case "$opt" in
	h)
	    show_help
	    exit 0
	    ;;
	c)
	    colindexfile=$OPTARG
	    ;;
	r)
	    rowindexfile=$OPTARG
	    ;;
	d)
	    countsfile=$OPTARG
	    ;;
	o)
	    outputfile=$OPTARG
	    ;;
	t)
	    tmpdir=$OPTARG
	    ;;
    esac
done

shift $((OPTIND-1))

## DEBUG check args
# echo colindex $colindexfile
# echo rowindex $rowindexfile
# echo counts $countsfile
# echo outputfile $outputfile

## Check the input
if [ ! -f $colindexfile ]; then
    echo "Column index file does not exist!";
    exit;
fi

if [ ! -f $rowindexfile ]; then
    echo "Row index file does not exist!";
    exit;
fi

if [ ! -f $countsfile ]; then
    echo "Counts file does not exist!";
    exit;
fi

if [ -f $outputfile ]; then
    echo "The output file already exists!";
    exit;
fi

## Make tmp directory
if [ -z "$tmpdir" ]; then
    tmpdir=`mktemp -d`
fi

if [ ! -d "$tmpdir" ]; then
    echo "Error: Temporary directory '$tmpdir' is not a directory"
    exit;
fi

## Convert the npz to text
echo "Converting npz to text..."
npz2txt.py --col-index $colindexfile \
	     --row-index $rowindexfile \
	     --counts $countsfile \
	     --output-dir $tmpdir

## Convert the text to rds
echo "Converting text to rds..."
sparseTxt2Rds.R --input-dir $tmpdir --output-file $outputfile

