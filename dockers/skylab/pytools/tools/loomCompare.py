#!/usr/bin/env python3

import sys
import argparse
import loompy
import numpy as np
import pandas as pd


def main():
    description = """This script compares two loom files and checks that they contain identical data up to a constant"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--truth-loom", dest="truth_loom_path", required=True, help="Path to truth loom file", type=str)
    parser.add_argument("--check-loom", dest="check_loom_path", required=True, help="Path to loom file to check", type=str)
    parser.add_argument("--delta-cutoff", dest="delta_cutoff", required=True, help="Max delta value allowed", type=int)
    args = parser.parse_args()
    truth_loom = loompy.connect(args.truth_loom_path)
    check_loom = loompy.connect(args.check_loom_path)

    truth_loom_array = pd.DataFrame(data=truth_loom[:, :], index=truth_loom.row_attrs['gene_names'], columns=truth_loom.col_attrs['cell_names'])
    check_loom_array = pd.DataFrame(data=check_loom[:,:], index=check_loom.row_attrs['gene_names'], columns=check_loom.col_attrs['cell_names'])
    check_loom_array = check_loom_array[truth_loom_array.columns]

    delta = (check_loom_array - truth_loom_array).abs().sum().sum()

    if delta < args.delta_cutoff:
        print(f"Matrices are identical: delta: {delta} delta_cutoff: {args.delta_cutoff}")
        sys.exit(0)
    else:
        print(f"Matrices are NOT identical: delta: {delta} delta_cutoff: {args.delta_cutoff}")
        sys.exit(1)


if __name__ == "__main__":
    main()

