#!/usr/bin/env python3

import os
import argparse
import loompy

def main():
    description = """Merge the outputs of multiple SS2 pipeline runs into a single Loom file"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-loom-files',
                        dest='input_loom_files',
                        nargs="+",
                        required=True,
                        help="Path to input loom directory in DirectoryStore format")
    parser.add_argument('--output-loom-file',
                        dest='output_loom_file',
                        required=True,
                        help="Path to output loom file")
    parser.add_argument('--plate-sample-id',
                        dest='plate_sample_id',
                        required=True,
                        help="Plate sample id for output loom")
    args = parser.parse_args()

    # The list of Loom files that we need to merge
    
    loom_file_list = args.input_loom_files
    
    attrDict = dict()
    attrDict['sample_id'] = args.plate_sample_id
    loompy.combine(loom_file_list,output_file=args.output_loom_file,file_attrs = attrDict)

if __name__ == '__main__':
    main()
