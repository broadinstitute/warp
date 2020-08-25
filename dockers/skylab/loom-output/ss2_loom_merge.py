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
    parser.add_argument('--batch_id',
                        dest='batch_id',
                        required=True,
                        help="Batch id for output loom")
    parser.add_argument('--batch_name',
                        dest='batch_name',
                        help='User provided plate id for output loom')
    parser.add_argument('--pipeline_version',
                        dest='pipeline_version',
                        required=True,
                        help='Multisample SS2 version')
    args = parser.parse_args()

    # The list of Loom files that we need to merge
    
    loom_file_list = args.input_loom_files
    
    attrDict = dict()
    attrDict['batch_id'] = args.batch_id
    attrDict['pipeline_version'] = args.pipeline_version
    if args.batch_name is not None:
        attrDict['batch_name'] = args.batch_name

    loompy.combine(loom_file_list,output_file=args.output_loom_file,file_attrs = attrDict)

if __name__ == '__main__':
    main()
