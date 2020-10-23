#!/usr/bin/env python3

import argparse
import loompy
import pandas as pd
import json

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
    # The list of Loom files that we need to merge
    
    df = pd.read_table(args.input_file)
    loom_file_list = df['loom_file_name'].values
    
    attrDict = dict()
    for k in df.columns:
        attrDict[k]=",".join(df[k].tolist())

    for i in range(df.shape[0]):
        jsonDict = dict()
        for j in range(df.shape[1]):
            jsonDict[str(df.columns[j])] = df[df.columns[j]][i]
        print(jsonDict)
        output_name = str(df['loom_file_name'][i]).replace('.loom', '') + ".json"
        with open(output_name, "w") as f:
            json.dump(jsonDict, f)
    loompy.combine(loom_file_list,output_file=args.output_loom_file, file_attrs=attrDict)

if __name__ == '__main__':
    main()

