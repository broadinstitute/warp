#!/usr/bin/env python3

import argparse
import loompy
import pandas as pd

def main():
    description = """Merge the outputs of multiple SS2 pipeline runs into a single Loom file"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-file',
                        dest='input_file',
                        required=True,
                        help="Path to input metadata tsv")
    parser.add_argument('--output-loom-file',
                        dest='output_loom_file',
                        required=True,
                        help="Path to output loom file")
    args = parser.parse_args()
    print(args)
    # The list of Loom files that we need to merge

    df = pd.read_table(args.input_file)
    print(df.shape)
    loom_file_list = df['loom_file_name'].values
    attrDict = {}
    for k in df.columns:
        attrDict[k]=df.
    print(loom_file_list)


    loompy.combine(loom_file_list,output_file=args.output_loom_file)

if __name__ == '__main__':
    main()