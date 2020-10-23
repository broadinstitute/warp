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

    keys = df.columns[1:]
    print(keys)

    for row_num in range(df.shape[0]):
        print(row_num)
        f = loom_file_list[row_num]
        print(f)
        ds = loompy.connect(f)

        for col_num in range(len(keys)):
            ds.ca[keys[col_num]] = [df.iloc[row_num][col_num + 1]]
        ds.close()

#    for i in range(df.shape[0]):
        jsonDict = dict()
        for j in range(df.shape[1]):
            column_name = df.columns[j]
            jsonDict[str(column_name)] = str(df[column_name][row_num])
        print(jsonDict)
        output_name = str(df['loom_file_name'][row_num]).replace('.loom', '') + ".json"
        with open(output_name, "w") as f:
            json.dump(jsonDict, f)
    loompy.combine(loom_file_list, output_file=args.output_loom_file)

if __name__ == '__main__':
    main()

