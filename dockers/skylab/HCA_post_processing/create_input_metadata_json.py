#!/usr/bin/env python3

import json
import argparse


def main():
    description = """Collects input metadata from individual analysis file jsons """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-json-files',
                        dest='input_files',
                        nargs="+",
                        required=True,
                        help="List of son files")
    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help="Name of output file")

    args = parser.parse_args()

    analysis_files = args.input_files

    inputs = {"inputs": []}

    for analysis_file in analysis_files:
        with open(analysis_file, "r") as f:
            analysis_metadata = json.load(f)
        if analysis_metadata["file_core"]["file_name"].endswith(".loom"):
            input_uuid = analysis_metadata["provenance"]["document_id"]
            inputs["inputs"].append({"input_id": input_uuid, "input_type": "analysis_file"})

    with open(args.output, "w") as f:
        json.dump(inputs, f)


if __name__ == '__main__':
    main()

