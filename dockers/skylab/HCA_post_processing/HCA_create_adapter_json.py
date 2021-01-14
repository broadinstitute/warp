#!/usr/bin/env python3

import argparse
import json
import uuid
import re
import os
import subprocess

NAMESPACE = uuid.UUID('c6591d1d-27bc-4c94-bd54-1b51f8a2456c')


def get_uuid5(value_to_hash):
    return str(uuid.uuid5(NAMESPACE, value_to_hash))


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
    description = """Creates json files needed for HCA DCP2 MVP"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--project-loom-file',
                        dest='project_loom_file',
                        required=True,
                        help="Path to project loom file")
    parser.add_argument('--crc32c',
                        dest='crc32c',
                        required=True,
                        help="crc32c of the loom file")
    parser.add_argument('--version-timestamp',
                        dest='version_timestamp',
                        required=True,
                        help="A version for the output files in the form of a timestamp")
    parser.add_argument('--project-id',
                        dest='project_id',
                        required=True,
                        help="project id of the loom file")
    parser.add_argument('--project-stratum-string',
                        dest='project_stratum_string',
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
    parser.add_argument('--loom-timestamp',
                        dest='loom_timestamp',
                        required=True,
                        help="The timestamp for the stratified project matrix loom file")
    parser.add_argument('--pipeline-version',
                        dest='pipeline_version',
                        required=True,
                        help="The version of the pipeline used to create the stratified project matrix")

    args = parser.parse_args()

    project_loom_file = args.project_loom_file
    crc32c = args.crc32c
    file_version = args.version_timestamp
    loom_timestamp = args.loom_timestamp
    loom_version = loom_timestamp.replace('Z', '.000000Z')
    project_id = args.project_id
    project_stratum_string = args.project_stratum_string
    sha256 = args.sha256
    size = int(args.size)
    staging_bucket = args.staging_bucket
    pipeline_version = args.pipeline_version
    with open(args.inputs_json, "r") as i:
        inputs_dict = json.load(i)  # this should be a list of dictionaries
        inputs = inputs_dict['inputs']

    analysis_type = "run"
    if "cacheCopy" in str(project_loom_file):
        analysis_type = "copy-forward"

    # Generate additional data from args
    file_name = os.path.basename(project_loom_file)
    process_id = get_analysis_workflow_id(project_loom_file)

    # Create UUIDs
    links_id = get_uuid5(project_stratum_string)  # v5 UUID of project id and the values the data are stratified by
    matrix_entity_id = get_uuid5(str(links_id + "analysis_file" + "loom"))   # v5 UUID of the links_id
    matrix_file_id = get_uuid5(matrix_entity_id)  # v5 UUID of the matrix_entity_id

    analysis_file_dict = {
                           "describedBy": "https://schema.humancellatlas.org/type/file/6.2.0/analysis_file",
                           "file_core": {
                             "file_name": file_name,
                             "format": "loom",
                             "content_description": [{
                               "text": "DCP/2-generated matrix",
                               "ontology": "data:3917",
                               "ontology_label": "Count Matrix"
                             }]
                           },
                           "provenance": {
                             "document_id": matrix_entity_id,
                             "submission_date": file_version,
                             "submitter_id": "e67aaabe-93ea-564a-aa66-31bc0857b707"
                           },
                           "schema_type": "file"
                         }

    analysis_process_dict = {
                              "describedBy": "https://schema.humancellatlas.org/type/process/analysis/12.0.0/analysis_process",
                              "schema_type": "process",
                              "process_core": {
                                "process_id": process_id
                              },
                              "type": {
                                "text": "analysis; merge matrices"
                              },
                              "reference_files": [],
                              "timestamp_start_utc": loom_version,  # string;
                                                                    # Initial start time of the full pipeline in UTC.
                                                                    # format: yyyy-mm-ddThh:mm:ssZ
                              "timestamp_stop_utc": loom_version,   # string;
                                                                    # Terminal stop time of the full pipeline in UTC.
                                                                    # format: yyyy-mm-ddThh:mm:ssZ
                              "tasks": [
                               # {
                               #  "task_name": "",   # string; Name of the task.
                               #                     # example: CollectDuplicationMetrics; RSEMExpression
                               #  "start_time": "",  # string; Date and time when the task started.
                               #                     # Enter the time in date-time format: yyyy-mm-ddThh:mm:ssZ
                               #  "stop_time": "",   # string; Date and time when the task finished.
                               #                     # Enter the time in date-time format: yyyy-mm-ddThh:mm:ssZ
                               #  "disk_size": "",   # string; Name of the disk volume mounted to the VM for the task.
                               #                     # Indicate both disk type and disk size. example: local-disk 11 HDD
                               #  "docker_image": "",# string;
                               #                     # Name of docker image where the task is stored and executed.
                               #                     # quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10
                               #  "cpus": 0,         # integer; Number of CPUs used to run this task.
                               #  "memory": "",      # string; Amount of memory allocated for this task. example: 7.5 GB
                               #  "zone": ""         # string Name of the Google Cloud zone where the task was run.
                               #                     #example: us-central1-b; europe-north1-a
                               # }
                              ],
                              "inputs": [
                               # {
                               #   "parameter_name": "",  # string; Name of parameter. example: stranded; rsem_ref_index
                               #   "parameter_value": ""  # string; Path to file for or value of parameter.
                               #                          # example: NONE;
                               #                          # gs://hca-dcp-mint-test-data/../gencode_v27_primary.tar"
                               # }  # Input parameters used in the pipeline run.
                              ],
                              "analysis_run_type": analysis_type,
                              "provenance": {
                                "document_id": process_id,
                                "submission_date": file_version,
                              },
                            }

    analysis_protocol_dict = {
                               "describedBy": "https://schema.humancellatlas.org/type/protocol/analysis/9.1.0/analysis_protocol",
                               "schema_type": "protocol",
                               "protocol_core": {
                                 "protocol_id": pipeline_version
                               },
                               "computational_method": pipeline_version,  # string; A URI to a versioned workflow and
                                                                          # versioned execution environment in a
                                                                          # GA4GH-compliant repository.
                                                                          # example: SmartSeq2SingleCell; 10x
                               "type": {
                                 "text": "analysis; merge matrices"
                               }
                             }
    analysis_protocol_string = json.dumps(analysis_protocol_dict, sort_keys=True)
    analysis_protocol_entity_id = get_uuid5(analysis_protocol_string)
    analysis_protocol_dict['provenance'] = {
                                             'document_id': analysis_protocol_entity_id,
                                             'submission_date': file_version,
                                             'update_date': file_version
                                           }

    file_descriptor_dict = {
                             "crc32c": crc32c,
                             "content_type": "application/vnd.loom",
                             "describedBy": "https://schema.humancellatlas.org/system/2.0.0/file_descriptor",
                             "file_id": matrix_file_id,
                             "file_name": file_name,
                             "file_version": loom_version,
                             "schema_type": "file_descriptor",
                             "schema_version": "2.0.0",
                             "sha256": sha256,
                             "size": size,
                           }

    links_dict = {
                   "describedBy": "https://schema.humancellatlas.org/system/2.1.1/links",
                   "links": [
                     {
                       "inputs": inputs,
                       "link_type": "process_link",
                       "outputs": [
                         {
                           "output_id": matrix_entity_id,
                           "output_type": "analysis_file"
                         }
                       ],
                       "process_id": process_id,
                       "process_type": "analysis_process",
                       "protocols": [
                         {
                           "protocol_id": analysis_protocol_entity_id,
                           "protocol_type": "analysis_protocol"
                         }
                       ]
                     }
                   ],
                   "schema_type": "links",
                   "schema_version": "2.1.1"
                 }

    # filenames for staging directories
    analysis_file_basename = "{}_{}.json".format(matrix_entity_id, file_version)
    analysis_protocol_basename = "{}_{}.json".format(analysis_protocol_entity_id, file_version)
    analysis_process_basename = "{}_{}.json".format(process_id, file_version)
    links_basename = "{}_{}_{}.json".format(links_id, file_version, project_id)

    # files created in output directory for output
    analysis_file_json_file_name = "outputs/analysis_file_{}".format(analysis_file_basename)
    analysis_process_json_file_name = "outputs/analysis_process_{}".format(analysis_process_basename)
    analysis_protocol_json_file_name = "outputs/analysis_protocol_{}".format(analysis_protocol_basename)
    file_descriptor_json_file_name = "outputs/file_descriptor_{}".format(analysis_file_basename)
    links_json_file_name = "outputs/links_{}".format(links_basename)

    with open(analysis_file_json_file_name, "w") as f:
        json.dump(analysis_file_dict, f, sort_keys=True, indent=2)

    with open(analysis_process_json_file_name, "w") as f:
        json.dump(analysis_process_dict, f, sort_keys=True, indent=2)

    with open(analysis_protocol_json_file_name, "w") as f:
        json.dump(analysis_protocol_dict, f, sort_keys=True, indent=2)

    with open(file_descriptor_json_file_name, "w") as f:
        json.dump(file_descriptor_dict, f, sort_keys=True, indent=2)

    with open(links_json_file_name, "w") as f:
        json.dump(links_dict, f, sort_keys=True, indent=2)

    # Copy json files into the staging bucket
    subprocess.run('gsutil cp {0} {1}data/{2}'.format(project_loom_file, staging_bucket, file_name), shell=True)
    subprocess.run('gsutil cp {0} {1}metadata/analysis_file/{2}'.format(analysis_file_json_file_name,
                                                                        staging_bucket,
                                                                        analysis_file_basename), shell=True)
    subprocess.run('gsutil cp {0} {1}metadata/analysis_process/{2}'.format(analysis_process_json_file_name,
                                                                           staging_bucket,
                                                                           analysis_process_basename), shell=True)
    subprocess.run('gsutil cp {0} {1}metadata/analysis_protocol/{2}'.format(analysis_protocol_json_file_name,
                                                                            staging_bucket,
                                                                            analysis_protocol_basename), shell=True)
    subprocess.run('gsutil cp {0} {1}descriptors/analysis_file/{2}'.format(file_descriptor_json_file_name,
                                                                           staging_bucket,
                                                                           analysis_file_basename), shell=True)
    subprocess.run('gsutil cp {0} {1}links/{2}'.format(links_json_file_name,
                                                       staging_bucket,
                                                       links_basename), shell=True)


if __name__ == '__main__':
    main()

