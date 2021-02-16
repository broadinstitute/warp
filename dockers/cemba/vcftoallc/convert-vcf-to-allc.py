#!/usr/bin/env python3
import argparse
import gzip
import re

def main():
    """ This script converts a VCF to ALLC file format
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
      "--input_vcf", "-i", dest="input_vcf", default=None, required=True, help="input VCF"
    )
    parser.add_argument(
      "--output_allc", "-o", dest="output_allc", default=None, help="output ALLC"
    )
    args = parser.parse_args()

    writeBinary = False
    with gzip.open(args.input_vcf, "rt") if args.input_vcf.endswith(".gz"
    ) else open(args.input_vcf, "r") as vcf_file:
        with gzip.open(args.output_allc, "wb") if args.output_allc.endswith(".gz"
        ) else open(args.output_allc, "w") as allc_file:

            writeBinary = args.output_allc.endswith(".gz")

            for line in vcf_file:
                if line.startswith("#"):
                    continue
                lineSplits = line.strip().split("\t")

                testInfo = lineSplits[7].split(";")
                testContext = re.sub("REFERENCE_CONTEXT=", "", testInfo[2])
                testUnMethCov = int(re.sub("CONVERTED_BASE_COV=", "", testInfo[0]))
                testMethCov = int(re.sub("UNCONVERTED_BASE_COV=", "", testInfo[3]))
                testTotCov = testMethCov + testUnMethCov

                ALLC = [lineSplits[0], lineSplits[1]]

                if lineSplits[3] == "G":
                    ALLC.append("-")
                else:
                    ALLC.append("+")
                ALLC.append(testContext)
                ALLC.append(str(testMethCov))
                ALLC.append(str(testTotCov))
                ALLC.append("1")

                if writeBinary:
                    allc_file.write(("\t".join(ALLC) + "\n").encode())
                else:
                    allc_file.write("\t".join(ALLC) + "\n")

if __name__ == "__main__":
    main()
