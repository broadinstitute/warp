#!/usr/bin/env python3

import argparse
import json
import uuid

NAMESPACE = uuid.UUID('c6591d1d-27bc-4c94-bd54-1b51f8a2456c')

def get_uuid5(sha256):
    return str(uuid.uuid5(NAMESPACE, sha256))

def main():
    description = """Add metadata into a Loom file as column attributes"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-loom-file',
                        dest='input_loom_file',
                        required=True,
                        help="Path to project loom file")
    parser.add_argument('--size',
                        dest='size',
                        required=True,
                        help="Size of the loom file in bytes")
    parser.add_argument('--sha256',
                        dest='sha256',
                        required=True,
                        help="sha256 of the loom file")
    parser.add_argument('--crc32c',
                        dest='crc32c',
                        required=True,
                        help="crc32c of the loom file")
    parser.add_argument('--version',
                        dest='version',
                        required=True,
                        help="version of the loom file")

    args = parser.parse_args()


    size = args.size
    sha256 = args.sha256
    crc32c = args.crc32c
    version = args.version

    file_version # The version of the file given in date time format
    file_name # The object name of the data file relative to the staging area's `data/` directory
    submission_date # When project was first submitted to database. date-time format # TODO ??????
    content_type = input_type = output_type = format ???
    process_id # UUID of the process described by this link # cromwell id of the process that generated the combined matrix TODO how?
    process_type # The concrete type of the process described by this link - "example": "analysis" TODO ask in channel if there are expected values
    protocols # An array of protocols for this link

    # TODO need a function to generate this for each matrix used as input
    inputs # need to create a list of inputs each with input_type and input_id defined

    # is format the same as input_type and output_type (and content_type ?)


    matrix_file_uuid = get_uuid5(sha256)
    file_uuid = get_uuid5(matrix_file_uuid)

    # need an analysis file instead of the supplementary file
    analysis_file_dict = {
                          "provenance": {
                                         "document_id": matrix_file_uuid,
                                         "submission_date": submission_date
                                        },
                          "describedBy": "https://schema.humancellatlas.org/type/file/6.2.0/analysis_file",
                          "schema_type": "file",
                          "file_core": {
                                        "file_name": file_name, # "The name of the file. Include the file extension in the file name."
                                        "format": "loom" # Indicate the full file extension including compression extensions.
                                       }
                         }

    file_descriptor_dict = {"describedBy": "https://schema.humancellatlas.org/system/2.0.0/file_descriptor",
                            "schema_type": "file_descriptor",
                            "content_type": content_type,
                            "size": size,
                            "sha256": sha256,
                            "crc32c": crc32c,
                            "file_id": file_uuid,
                            "file_version": file_version,
                            "file_name": file_name
                            }

    links_dict = {"describedBy": "https://schema.humancellatlas.org/system/2.1.1/links",
                  "schema_type": "links",
                  "links": [ { "process_type": process_type,
                               "process_id": process_id,
                               "inputs": [ {
                                            "input_type": "loom",
                                            "input_id": "UUID"
                                            }
                                         ],
                               "outputs": [ {
                                            "output_type": "loom",
                                            "output_id": matrix_file_uuid
                                            }
                                         ],
                               "protocols": protocols,
                               "link_type": "process_link"
                              }
                            ]
                  }

    # file_version is a timestamp and will be the same for all instances below

    analysis_file_json_file_name = "metadata/{}/{}_{}.json".format(entity_type, matrix_file_uuid, file_version)
    file_descriptor_json_file_name = "descriptors/{}/{}_{}.json".format(entity_type, matrix_file_uuid, file_version)
    links_json_file_name = "links/{}_{}_{}.json".format(matrix_file_uuid, file_version, project_id)

    # TODO move the matrix file to data/{file_name}

    with open(analysis_file_json_file_name, "w") as f:
        json.dump(analysis_file_dict, f)

    with open(file_descriptor_json_file_name, "w") as f:
        json.dump(file_descriptor_dict, f)

    with open(links_json_file_name, "w") as f:
        json.dump(links_dict, f)

if __name__ == '__main__':
    main()

