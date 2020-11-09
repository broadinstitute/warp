#!/usr/bin/env python3

import json
import argparse


def main():
    description = """Add metadata into a Loom file as column attributes"""
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

    inputs = {"inputs": []}

    for input_file in args.input_files:
        analysis_metadata = json.loads(input_file)
        if analysis_metadata["file_core"]["file_name"].endswith("*.loom"):
            input_uuid = analysis_metadata["provenance"]["document_id"]

        inputs["inputs"].append({"input_id": input_uuid, "input_type": "analysis_file"})

    with open(args.output, "w") as f:
        json.dump(inputs, f)


if __name__ == '__main__':
    main()

