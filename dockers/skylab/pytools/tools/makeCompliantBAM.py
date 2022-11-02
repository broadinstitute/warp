#!/usr/bin/env python

import pysam
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Make a BAM file with cellular barcodes in the "
                                                 + "read names compliant by moving them to the CB tag")
    parser.add_argument('--input-bam', dest="inputbam", help="input bam file")
    parser.add_argument('--output-bam', dest="outputbam", help="output bam file")

    args = parser.parse_args()

    def checkArgs(args):
        if args.inputbam is None:
            sys.exit("Input BAM is not defined")
        if args.outputbam is None:
            sys.exit("Output BAM is not defined")
        if not os.path.isfile(args.inputbam):
            sys.exit("Input BAM is not a file")

    checkArgs(args);

    inputbamfilename = "input.bam"
    bamfile = pysam.AlignmentFile(args.inputbam, 'rb')

    outputfilename = "output.bam"
    outbam = pysam.AlignmentFile(args.outputbam, 'wb', template=bamfile)

    counter = 0

    for read in bamfile:
        counter += 1
        if (counter % 100000 == 0):
            print(counter)
        qname = str(read.qname)
        i = qname.find(':')
        cb, qn = qname[:i], qname[i+1:]
        read.qname = qn
        read.set_tag('CB',cb,'Z')
        outbam.write(read)

    outbam.close()

    bamfile.close()

if __name__ == '__main__':
    main()

