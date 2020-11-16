#!/usr/bin/env python3

import argparse
import json
import uuid
import re
import os
import subprocess

NAMESPACE = uuid.UUID('c6591d1d-27bc-4c94-bd54-1b51f8a2456c')

def get_uuid5(sha256):
    return str(uuid.uuid5(NAMESPACE, sha256))

def get_analysis_workflow_id(analysis_output_path):
    """Parse the analysis workflow id from one of its output paths, and write the id to a file so that it is available
    outside of the get_analysis task.
    Args:
        analysis_output_path (str): path to workflow output file.
    Returns:
        workflow_id (str): string giving Cromwell UUID of the workflow.
    """
    # Get the last match for UUID prior to the file name (in case the file is
    # named with a UUID) to ensure it is the subworkflow id
    url = analysis_output_path.rsplit('/', 1)[0]
    uuid_regex = r"([a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12})"
    workflow_id = re.findall(uuid_regex, url)[-1]
    print('Got analysis workflow UUID: {0}'.format(workflow_id))
    return workflow_id


def main():
    description = """Add metadata into a Loom file as column attributes"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--project-loom-file',
                        dest='project_loom_file',
                        required=True,
                        help="Path to project loom file")
    parser.add_argument('--crc32c',
                        dest='crc32c',
                        required=True,
                        help="crc32c of the loom file")
    parser.add_argument('--file_version',
                        dest='file_version',
                        required=True,
                        help="timepstamp of the loom file")
    parser.add_argument('--project-id',
                        dest='project_id',
                        required=True,
                        help="project id of the loom file")
    parser.add_argument('--sha256',
                        dest='sha256',
                        required=True,
                        help="sha256 of the loom file")
    parser.add_argument('--size',
                        dest='size',
                        required=True,
                        help="Size of the loom file in bytes")
    parser.add_argument('--staging-bucket',
                        dest='staging_bucket',
                        help="Path to staging bucket")
    parser.add_argument('--input-metadata-json',
                        dest='inputs_json',
                        required=True,
                        help="Json file with inputs metadata")
    parser.add_argument('--protocol-metadata-json',
                        dest='protocols_json',
                        required=True,
                        help="Json file with protocols metadata")

    args = parser.parse_args()

    project_loom_file = args.project_loom_file
    crc32c = args.crc32c
    file_version = args.file_version
    project_id = args.project_id
    sha256 = args.sha256
    size = args.size
    staging_bucket = args.staging_bucket
    with open(args.inputs_json, "r") as i:
        inputs_dict = json.load(i)  # this should be a list of dictionaries
        inputs = inputs_dict['inputs']
    with open(args.protocols_json, "r") as p:
        protocols_dict = json.load(p)  # this should be a list of dictionaries
        protocols = protocols_dict['protocols']

    # Generate additional data from agrs
    file_name = os.path.basename(project_loom_file)
    process_id = get_analysis_workflow_id(project_loom_file)

    matrix_file_uuid = get_uuid5(sha256)
    file_uuid = get_uuid5(matrix_file_uuid)

    analysis_file_dict = {
                          "provenance": {
                                         "document_id": matrix_file_uuid,
                                         "submission_date": file_version
                                        },
                          "describedBy": "https://schema.humancellatlas.org/type/file/6.2.0/analysis_file",
                          "schema_type": "file",
                          "file_core": {
                                        "file_name": file_name,
                                        "format": "loom"
                                       }
                         }

    file_descriptor_dict = {"describedBy": "https://schema.humancellatlas.org/system/2.0.0/file_descriptor",
                            "schema_type": "file_descriptor",
                            "content_type": "application/vnd.loom",
                            "size": size,
                            "sha256": sha256,
                            "crc32c": crc32c,
                            "file_id": file_uuid,
                            "file_version": file_version,
                            "file_name": file_name
                            }

    links_dict = {"describedBy": "https://schema.humancellatlas.org/system/2.1.1/links",
                  "schema_type": "links",
                  "links": [{"process_type": "analysis",
                             "process_id": process_id,
                             "inputs": inputs,
                             "outputs": [{
                                         "output_id": matrix_file_uuid,
                                         "output_type": "analysis_file"
                                         }],
                             "protocols": protocols,
                             "link_type": "process_link"
                             }]
                  }

    # filenames for staging dierctories
    file_basename = "{}_{}".format(matrix_file_uuid, file_version)
    links_basename = "{}_{}_{}".format(matrix_file_uuid, file_version, project_id)

    # files created in output directory for output
    analysis_file_json_file_name = "outputs/analysis_file_{}.josn".format(file_basename)
    file_descriptor_json_file_name = "outputs/file_descriptor_{}.json".format(file_basename)
    links_json_file_name = "outputs/links_{}.json".format(links_basename)

    with open(analysis_file_json_file_name, "w") as f:
        json.dump(analysis_file_dict, f)

    with open(file_descriptor_json_file_name, "w") as f:
        json.dump(file_descriptor_dict, f)

    with open(links_json_file_name, "w") as f:
        json.dump(links_dict, f)

    # Copy json files into the staging bucket
    subprocess.run('gsutil cp {0} {1}/data/{0}'.format(project_loom_file, staging_bucket), shell=True)
    subprocess.run('gsutil cp {0} {1}/metadata/analysis_file/{2}.json'.format(analysis_file_json_file_name,
                                                                         staging_bucket,
                                                                         file_basename), shell=True)
    subprocess.run('gsutil cp {0} {1}/descriptors/file_descriptor/{2}.json'.format(file_descriptor_json_file_name,
                                                                              staging_bucket,
                                                                              file_basename), shell=True)
    subprocess.run('gsutil cp {0} {1}/links/{2}.json'.format(links_json_file_name,
                                                        staging_bucket,
                                                        links_basename), shell=True)


if __name__ == '__main__':
    main()

