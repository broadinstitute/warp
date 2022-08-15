#!/usr/bin/env python

import argparse
import h5py
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Extract the data in a snap file as csv files")
    parser.add_argument('--input-file', dest='inputfile', help='snap input file')
    parser.add_argument('--output-prefix',dest='output_prefix', help='output file prefix',default='')
    args = parser.parse_args()

    #inputfile = 'atac_v1_adult_brain_fresh_5k.snap'
    #hd5file = h5py.File('output-big.snap')

    hd5file = h5py.File(args.inputfile)

    # HD = hd5file['HD']
    # if not HD['MG'][()] == 'SNAP':
    #    raise('The input is not a SNAP file')

    # ###########################
    # Export the fragment matrix (FM)
    # This has the following columns: fragChrom, fragStart, fragLen, barcodePos, barcodeLen
    FM = hd5file['FM']
    # Export barcodes
    barcode_DF = pd.DataFrame(
        {
            'barcodeLen': FM['barcodeLen'][:],
            'barcodePos': FM['barcodePos'][:]
         }
    )
    barcode_DF.to_csv(args.output_prefix + 'barcodes.csv')

    # Export fragments
    fragments_DF = pd.DataFrame(
        {
            'fragChrom': FM['fragChrom'][:],
            'fragLen': FM['fragLen'][:],
            'fragStart': FM['fragStart'][:]
        }
    )
    fragments_DF.to_csv(args.output_prefix + 'fragments.csv')

    ###################################
    # Export the cell x bin accessibility matrices (AM)
    AM=hd5file['AM']

    # Snap files can contain multiple different summarizations into
    # bins of different sizes. We extract them all here
    nBinSize = AM['nBinSize'][()]
    bins = AM['binSizeList'][:]

    for binname in bins:
        binname = str(binname)
        binContainer = AM[binname]
        ## binCoordinates
        binCoordinates_DF = pd.DataFrame(
            {
                'binChrom': binContainer['binChrom'],
                'binStart': binContainer['binStart']
            }
        )
        binCoordinates_DF.to_csv(args.output_prefix + 'binCoordinates_' + str(binname) + '.csv')
        ## binCounts
        binCounts_DF = pd.DataFrame(
            {
                'idx': binContainer['idx'],
                'idy': binContainer['idy'],
                'count': binContainer['count']
            }
        )
        binCounts_DF.to_csv(args.output_prefix + 'binCounts_' + str(binname) + '.csv')

    ###################################
    # Extract the barcode session
    BD=hd5file['BD']
    barcodesSection_DF = pd.DataFrame(
        {
            'name': BD['name'],
            'TN': BD['TN'],
            'UM': BD['UM'],
            'SE': BD['SE'],
            'SA': BD['SA'],
            'PE': BD['PE'],
            'PP': BD['PP'],
            'PL': BD['PL'],
            'US': BD['US'],
            'UQ': BD['UQ'],
            'CM': BD['CM']
        }
    )
    barcodesSection_DF.to_csv(args.output_prefix + 'barcodesSection.csv')


if __name__ == '__main__':
    main()
