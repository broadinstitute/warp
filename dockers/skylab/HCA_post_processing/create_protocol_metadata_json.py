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

    all_protocols = {"protocols": []}
    ids = set([])

    for input_file in args.input_files:
        links_metadata = json.loads(input_file)
        protocols = links_metadata["links"]["protocols"]
        for protocol in protocols:
            if protocol["protocol_id"] not in ids:
                ids.add(protocol["protocol_id"])
                all_protocols["protocols"].append(protocol)

    with open(args.output, "w") as f:
        json.dump(all_protocols, f)


if __name__ == '__main__':
    main()

