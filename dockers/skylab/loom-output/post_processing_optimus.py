#!/usr/bin/env python3

import argparse
import loompy
import pandas as pd
import json

def main():
    description = """Merge the outputs of multiple optimus pipeline runs into a single Loom file"""
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
    metadata_df = df.copy()
    metadata_df = metadata_df.drop(columns='loom_file_name')
    metadata_df.drop_duplicates(inplace=True)

    if metadata_df.shape[0] > 1:
        raise ValueError('Incompatible metadata')

    keys = df.columns[1:]
    print(keys)

    for row_num in range(df.shape[0]):
        print(row_num)
        f = loom_file_list[row_num]
        print(f)
        ds = loompy.connect(f)

        ds.ca['cell_names'] = ds.ca['cell_names'] + "-" + str(row_num)
        print(ds.ca['cell_names'])

        ds.close()
        # TODO: add some map that give the sample/library that corresponds to the barcode suffix


        jsonDict = dict()
        for j in range(df.shape[1]):
            column_name = df.columns[j]
            jsonDict[str(column_name)] = str(df[column_name][row_num])
        print(jsonDict)
        output_name = str(df['loom_file_name'][row_num]).replace('.loom', '') + ".json"
        with open(output_name, "w") as f:
            json.dump(jsonDict, f)
    attr_dict = {k: v for k in metadata_df.columns for v in metadata_df[k]}
    # TODO: add map from above to attr_dict
    print(attr_dict)


    loompy.combine(loom_file_list, output_file=args.output_loom_file, file_attrs=attr_dict)

if __name__ == '__main__':
    main()

