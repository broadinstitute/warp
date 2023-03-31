#! /usr/bin/env python

# Script to append R2 barcodes to ATAC read 1 and read 3 FASTQs
from argparse import ArgumentParser
from itertools import islice

def write_fastq(args):
    """ Appends barcodes from R2 to R1 and R3 ATAC fastq files. 

        Args:
            input_r1_fastq (str): input R1 ATAC fastq file.
            input_r3_fastq (str): input R3 ATAC fastq file.
            cb_fastq (str): input R2 ATAC fastq file which contains the cell barcodes.
            out_r1_fastq (str): output file name for R1 fastq file with appended barcodes.   
            out_r3_fastq (str): output file name for R3 fastq file with appended barcodes.  
    """
    w_r1 = open(args.out_r1_fastq,"a+")
    w_r3 = open(args.out_r3_fastq, "a+")
    with open(args.input_r1_fastq, 'r+') as R1, open(args.cb_fastq, 'r+') as R2, open(args.input_r3_fastq, 'r+') as R3:
        while True:
            lines_R1=list(islice(R1, 4))
            lines_R2=list(islice(R2, 4))
            lines_R3=list(islice(R3, 4))
            if not lines_R1:
                break
            if lines_R1[0].startswith("@") and lines_R2[0].startswith("@") and lines_R3[0].startswith("@"):
                if lines_R1[0].split(" ")[0]  == lines_R2[0].split(" ")[0] == lines_R3[0].split(" ")[0]:
                    barcode = lines_R2[1].strip()
                    new_header_r1 = "@" + barcode + ":" + lines_R1[0][1:].strip()
                    new_header_r3 = "@" + barcode + ":" + lines_R3[0][1:].strip()
                    w_r1.write(new_header_r1 + "\n" + lines_R1[1].strip() + "\n" + lines_R1[2].strip() + "\n" + lines_R1[3].strip() + "\n")
                    w_r3.write(new_header_r3 + "\n" + lines_R3[1].strip() + "\n" + lines_R3[2].strip() + "\n" + lines_R3[3].strip() + "\n")
    w_r1.close()
    w_r3.close()


def main():
    parser = ArgumentParser(description = 'Appending cell barcodes to ATAC FASTQs')       
    parser.add_argument('-r1', dest= "input_r1_fastq", required = True, metavar='fastq with input reads', help='Input R1 FASTQ')
    parser.add_argument('-r3', dest= "input_r3_fastq", required = True, metavar='fastq with input reads', help='Input R3 FASTQ')
    parser.add_argument('-cb', dest= "cb_fastq", required = True, metavar='fastq with barcodes', help='Input R2 FASTQ with barcodes')
    parser.add_argument('-out_r1', dest= "out_r1_fastq", required = True, metavar='Output R1 FASTQ', help='Output FASTQ')
    parser.add_argument('-out_r3', dest= "out_r3_fastq", required = True, metavar='Output R3 FASTQ', help='Output FASTQ')
    args = parser.parse_args()


    write_fastq(args)
    

if __name__ == '__main__':
    main()

