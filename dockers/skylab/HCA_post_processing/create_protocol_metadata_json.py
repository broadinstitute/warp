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
    links_json_files = args.input_files

    all_protocols = {"protocols": []}
    ids = set([])

    for links_file in links_json_files:
        with open(links_file, "r") as f:
            links_metadata = json.load(f)
        protocols = links_metadata["links"]["protocols"]
        for protocol in protocols:
            if protocol["protocol_id"] not in ids:
                ids.add(protocol["protocol_id"])
                all_protocols["protocols"].append(protocol)

    with open(args.output, "w") as f:
        json.dump(all_protocols, f)


if __name__ == '__main__':
    main()

