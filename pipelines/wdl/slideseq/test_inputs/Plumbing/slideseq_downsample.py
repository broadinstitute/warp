#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import numpy as np


def main():
    description = """This script downsamples data for the SlideSeq pipeline for use in testing."""
    parser = argparse.ArgumentParser(description=description)

    # Input file arguments
    parser.add_argument("--xy_coords", dest="xy_coords", required=True, help="Path to bead xy-coordinates file",
                        type=str)
    parser.add_argument("--reads_umi", dest="reads_umi", required=True, help="Path to reads per UMI file", type=str)
    parser.add_argument("--reads_cbc", dest="reads_cbc", required=True, help="Path to reads per cell barcode", type=str)
    parser.add_argument("--dexpr_sum", dest="dexpr_sum", required=True, help="Path to digital expression summary file",
                        type=str)
    parser.add_argument("--raw_cbcs", dest="raw_cbcs", required=True, help="Path to raw cell barcodes file", type=str)

    # Cropping and sampling arguments
    # Test sample cropped from X: 2500 to 4000, Y: 1000 to 4200, sampling 10 cells per chunk in a 30 by 30 chunk region
    parser.add_argument("--xcoor_min", dest="xcoor_min", required=True, help="Minimum x-coordinate for sample cropping",
                        type=int)
    parser.add_argument("--xcoor_max", dest="xcoor_max", required=True, help="Maximum x-coordinate for sample cropping",
                        type=int)
    parser.add_argument("--ycoor_min", dest="ycoor_min", required=True, help="Minimum y-coordinate for sample cropping",
                        type=int)
    parser.add_argument("--ycoor_max", dest="ycoor_max", required=True, help="Maximum y-coordinate for sample cropping",
                        type=int)
    parser.add_argument("--cells_per_chunk", dest="cells_per_chunk", required=True, help="Cells per chunk to sample",
                        type=int)
    parser.add_argument("--nx_slices", dest="nx_slices", required=True,
                        help="Number of slices to sample along the x-axis", type=int)
    parser.add_argument("--ny_slices", dest="ny_slices", required=True,
                        help="Number of slices to sample along the y-axis", type=int)

    # Parse and store arguments
    args = parser.parse_args()

    # Create dataframes for input files
    xy_coord = pd.read_csv(args.xy_coords, sep="\t", header=None)
    reads_umi = pd.read_csv(args.reads_umi, sep="\t")
    nReads_CB = pd.read_csv(args.reads_cbc, sep="\t")
    expr_summary = pd.read_csv(args.dexpr_sum, sep="\t", skiprows=6)
    raw_barcodes = pd.read_csv(args.raw_cbcs, sep="\t", header=None, skiprows=1)

    # Clean up and merge dataframes
    xy_coord.columns = ["CELL_BARCODE", "x_coor", "y_coor"]
    expr_coor = expr_summary.merge(xy_coord)
    expr_coor['log_NUM_TRANSCRIPTS'] = np.log10(1 + expr_coor.NUM_TRANSCRIPTS)
    raw_barcodes.columns = ["CELL_BARCODE_RAW", "x_coor", "y_coor"]

    CB_counts = nReads_CB[nReads_CB.columns[0:2]]
    CB_counts.columns = ["nReads", "CELL_BARCODE"]
    CB_counts["log_nReads"] = np.log10(1 + CB_counts.nReads)

    raw_counts = xy_coord.merge(CB_counts)

    # Crop an interesting region from the sample
    filtered_counts = raw_counts.loc[(raw_counts.x_coor > args.xcoor_min) &
                                     (raw_counts.x_coor < args.xcoor_max) &
                                     (raw_counts.y_coor > args.ycoor_min) &
                                     (raw_counts.y_coor < args.ycoor_max)]

    # Downsample the cropped region by sampling chunks and write to a TSV file
    sample_set = filtered_counts.sample(1)
    xsteps = (args.xcoor_max - args.xcoor_min) / args.nx_slices
    ysteps = (args.ycoor_max - args.ycoor_min) / args.ny_slices

    for i in range(args.nx_slices):
        xlow = args.xcoor_min + xsteps * i
        xhigh = args.xcoor_min + xsteps * (i + 1)
        for j in range(args.ny_slices):
            ylow = args.ycoor_min + ysteps * j
            yhigh = args.ycoor_min + ysteps * (j + 1)
            subset_counts = filtered_counts.loc[(filtered_counts.x_coor > xlow) & (filtered_counts.x_coor < xhigh) &
                                                (filtered_counts.y_coor > ylow) & (filtered_counts.y_coor < yhigh)]
            sample_set = sample_set.append(subset_counts.sample(args.cells_per_chunk))

    sample_set.to_csv("sample.tsv", index=False, sep="\t")

    # Add raw barcodes to sample set and write to a TSV file
    barcodes_raw = []
    for x in sample_set.CELL_BARCODE:
        barcodes_raw.append(re.sub(r'-\d+$', '', x))
    sample_set['CELL_BARCODE_RAW'] = barcodes_raw
    matched_sample = sample_set.merge(raw_barcodes)

    matched_sample.to_csv("sample_with_barcodes.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()