import argparse
import os
import sys

import firecloud.api as fapi
from firecloud import fiss

import json_templates
import parse_terra_tsv as parse_terra

if "WORKSPACE_NAMESPACE" in os.environ:
    billing_project = os.environ['WORKSPACE_NAMESPACE']

template_workspace_name = "HCA_Smart-seq2_template"


def create_newworkspace(billing_project, template_workspace_name, new_workspace_name):
    ent_types = fiss.fapi.list_entity_types(billing_project, template_workspace_name)

    res = fiss.fapi.clone_workspace(billing_project,
                                    template_workspace_name,
                                    billing_project,
                                    new_workspace_name
                                    )

    return (res)


def get_data_from_workspace(billing_project, workspace_name, participant_table_name, output_table_name):
    r = fapi.get_entities_tsv(billing_project, workspace_name, participant_table_name)
    with open(output_table_name, 'w') as fout:
        fout.write(r.content.decode())

    if r.status_code != 200:
        print("ERROR :" + updated_workflow.content)
        sys.exit(1)
    else:
        print("Downloaded the table successfully")


def upload_tables(input_file, billing_project, workspace_name):
    with open(input_file) as tsvf:
        headerline = tsvf.readline().strip()
        entity_data = [l.rstrip('\n') for l in tsvf]

    model = 'flexible'

    batch_load(billing_project, workspace_name, headerline, entity_data, 500, model)


def upload_workflow_json(billing_project, workspace_name, namespace, workflow_name, json_templates):
    work_space_config = fapi.get_workspace_config(billing_project, workspace_name, namespace, workflow_name)
    work_space_json = work_space_config.json()
    work_space_json['inputs'] = json_templates.optimus_inputs
    work_space_json['outputs'] = json_templates.optimus_outputs

    updated_workflow = fapi.update_workspace_config(billing_project, workspace_name, namespace, workflow_name,
                                                    work_space_json)

    if updated_workflow.status_code != 200:
        print("ERROR :" + updated_workflow.content)
        sys.exit(1)
    else:
        print("updated successfully")


def create_job_submission(billing_project, workspace_name, workflow_name, workflow_repo):
    response = fapi.get_entities_with_type(billing_project, workspace_name)
    entities = response.json()

    for ent in entities:
        ent_name = ent['name']
        ent_type = ent['entityType']
        ent_attrs = ent['attributes']

    submisson_response = fapi.create_submission(billing_project, workspace_name,
                                                workflow_repo, workflow_name,
                                                entity=args.entity_id, etype="participant_lane_set",
                                                expression=None, use_callcache=True)

    if submisson_response.status_code != 201:
        print(submisson_response.content)
        sys.exit(1)
    else:
        print("Successfully Created Submisson")
        with open('response.txt', 'w') as fout:
            fout.write(submisson_response.json()['submissionId'] + '\n')


def valid_headerline(l, model='firecloud'):
    """return true if the given string is a valid loadfile header"""

    if not l:
        return False
    headers = l.split('\t')
    first_col = headers[0]

    tsplit = first_col.split(':')
    if len(tsplit) != 2:
        return False

    if tsplit[0] in ('entity', 'update'):
        if model == 'flexible':
            return tsplit[1].endswith('_id')
        else:
            return tsplit[1] in ('participant_id', 'participant_set_id',
                                 'participant_lane_id',
                                 'sample_id', 'sample_set_id',
                                 'pair_id', 'pair_set_id')
    elif tsplit[0] == 'membership':
        if len(headers) < 2:
            return False
        # membership:sample_set_id   sample_id, e.g.
        return tsplit[1].replace('set_', '') == headers[1]
    else:
        return False


def batch_load(project, workspace, headerline, entity_data, chunk_size=500,
               model='flexible'):
    """ Submit a large number of entity updates in batches of chunk_size """

    # Parse the entity type from the first cell, e.g. "entity:sample_id"
    # First check that the header is valid
    # if not valid_headerline(headerline, model):
    #    print("Invalid loadfile header:\n" + headerline)
    #    return 1

    update_type = "membership" if headerline.startswith("membership") else "entity"
    etype = headerline.split('\t')[0].split(':')[1].replace("_id", "")

    # Split entity_data into chunks
    total = int(len(entity_data) / chunk_size) + 1
    batch = 0
    for i in range(0, len(entity_data), chunk_size):
        batch += 1
        print("Updating {0} {1}s {2}-{3}, batch {4}/{5}".format(
            etype, update_type, i + 1, min(i + chunk_size, len(entity_data)),
            batch, total))
        this_data = headerline + '\n' + '\n'.join(entity_data[i:i + chunk_size])

        # Now push the entity data to firecloud
        r = fapi.upload_entities(project, workspace, this_data, model)
        fapi._check_response_code(r, 200)

    return 0


def main():
    if len(sys.argv) < 2:
        return

    global billing_project
    global template_workspace_name

    parser = argparse.ArgumentParser(prog="python " + sys.argv[0], add_help=False)
    subparser = parser.add_subparsers(dest="cmd")

    delete_workspace = subparser.add_parser('delete_workspace', help='delete workspace')
    delete_workspace.add_argument('--workspace-name', dest="workspace_name", help="name of the workspace")

    clone_workspace = subparser.add_parser('clone_workspace', help='clone from existing workspace')
    clone_workspace.add_argument('--source-work-space', dest='src_work_space', help="name of source workspace")
    clone_workspace.add_argument('--destination-work-space', dest='dest_work_space',
                                 help="name of destination workspace")

    get_data_info = subparser.add_parser('get_participant_table', help='get participant.tsv')
    get_data_info.add_argument('--workspace-name', dest="workspace_name", help="name of the workspace")
    get_data_info.add_argument('--participant-table-name', dest="participant_table_name", help="name of sample table")
    get_data_info.add_argument('--output-name', dest="output_table_name", required=False,
                               default="participant.tsv", help="name of output tsv"
                               )

    create_participant_lane = subparser.add_parser('create_participant_lane',
                                                   help='create participant_lane/lane_set_id tables')
    create_participant_lane.add_argument('--input-name', dest="input_participant_table_name",
                                         required=False,
                                         default="participant.tsv", help="input participant table  name")

    create_participant_lane.add_argument('--output-prefix', dest="output_prefix", required=False,
                                         help="name of output prefix for the lanes")

    upload_participant_lane = subparser.add_parser('upload_participant',
                                                   help='uploads the participant_lane_set, _lane_membership and _lane_entity files')
    upload_participant_lane.add_argument('--workspace-name', dest="workspace_name", help="name of the workspace")
    upload_participant_lane.add_argument('--input-prefix', dest="input_prefix", help="name of the input prefix")

    upload_workflow = subparser.add_parser('upload_workflow',
                                           help='uploads wdl to --workspace-name')
    upload_workflow.add_argument('--workspace-name', dest="workspace_name", help="name of the workspace")
    upload_workflow.add_argument('--method', dest="method", help="name of the input prefix")
    upload_workflow.add_argument('--wdl', dest="wdl", help="name of the input prefix")

    upload_config = subparser.add_parser('upload_config', help='upload config information')
    upload_config.add_argument('--workspace-name', dest="workspace_name", help="name of the workspace")
    upload_config.add_argument('--chemistry', dest="chemistry", choices=["V2", "V3"], help="chemistry")
    upload_config.add_argument('--counting-mode', dest="counting_mode", choices=["sc_rna", "sn_rna"],
                               help="counting mode: whether to count intronic alignments")
    upload_config.add_argument('--species', dest="species", choices=["human", "mouse"], help="species")

    submit_workflow = subparser.add_parser('submit_workflow', help='submit a workflow run')
    submit_workflow.add_argument('--workspace-name', dest="workspace_name", help="name of the workspace")
    submit_workflow.add_argument('--workflow-repo', dest="workflow_repo", help="workflow repo name")
    submit_workflow.add_argument('--workflow-name', dest="workflow_name", help="workflow name")
    submit_workflow.add_argument('--entity-id', dest="entity_id", help="entity id")

    get_status = subparser.add_parser('get_status', help='get status of a submission')
    get_status.add_argument('--workspace-name', dest="workspace_name", help="name of the workspace")
    get_status.add_argument('--submission-id', dest="submission_id", help="submission_id")

    # show help when no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # new_workspace_name = "DCP2_Optimus_template_KMK_v1"
    if args.cmd == 'delete_workspace':
        print("Delete existing workspace ", args.workspace_name)
        delete_status = fapi.delete_workspace(billing_project, args.workspace_name)

    elif args.cmd == 'clone_workspace':
        print("Cloning a new workspace from template", args.src_work_space)
        status = create_newworkspace(billing_project, args.src_work_space, args.dest_work_space)

    elif args.cmd == 'get_participant_table':
        print("Get information from workspace", args.workspace_name)
        r = fapi.get_entities_tsv(billing_project, args.workspace_name, args.participant_table_name)
        with open(args.output_table_name, 'w') as fout:
            fout.write(r.content.decode())

    elif args.cmd == 'create_participant_lane':
        parse_terra.create_output_files(args.input_participant_table_name, args.output_prefix)

    elif args.cmd == 'upload_participant':
        upload_tables(args.input_prefix + ".tsv", billing_project, args.workspace_name)
        upload_tables(args.input_prefix + "_membership.tsv", billing_project, args.workspace_name)
        upload_tables(args.input_prefix + "_entity.tsv", billing_project, args.workspace_name)
    elif args.cmd == 'upload_workflow':
        r = fapi.update_repository_method(args.workspace_name, args.method,
                                          "args.synopsis", args.wdl, "comment.txt",
                                          "args.comment")
        with open("response.txt", 'w') as fout:
            fout.write(r.content.decode())

    elif args.cmd == 'upload_config':

        work_space_config = fapi.get_workspace_config(billing_project, args.workspace_name, args.workspace_name,
                                                      "Optimus")
        work_space_json = work_space_config.json()
        work_space_json['inputs'] = json_templates.optimus_inputs
        work_space_json['outputs'] = json_templates.optimus_outputs

        if args.chemistry == "V2":
            work_space_json['inputs']['Optimus.chemistry'] = '\"tenX_v2\"'
            work_space_json['inputs']['Optimus.whitelist'] = 'workspace.whitelist_v2'
        if args.chemistry == "V3":
            work_space_json['inputs']['Optimus.chemistry'] = '\"tenX_v3\"'
            work_space_json['inputs']['Optimus.whitelist'] = 'workspace.whitelist_v3'

        if args.chemistry == "sn_rna":
            work_space_json['inputs']['Optimus.counting_mode'] = "\"sn_rna\""
        if args.chemistry == "sc_rna":
            work_space_json['inputs']['Optimus.counting_mode'] = "\"sc_rna\""

        if args.species == "human":
            work_space_json['inputs']['Optimus.annotations_gtf'] = 'workspace.human_annotations_gtf'
            work_space_json['inputs']['Optimus.ref_genome_fasta'] = 'workspace.human_ref_genome_fasta'
            work_space_json['inputs']['Optimus.tar_star_reference'] = 'workspace.human_tar_star_reference'
        if args.species == "mouse":
            work_space_json['inputs']['Optimus.annotations_gtf'] = 'workspace.mouse_annotations_gtf'
            work_space_json['inputs']['Optimus.ref_genome_fasta'] = 'workspace.mouse_ref_genome_fasta'
            work_space_json['inputs']['Optimus.tar_star_reference'] = 'workspace.mouse_tar_star_reference'

        updated_workflow = fapi.update_workspace_config(billing_project, args.workspace_name, args.workspace_name,
                                                        "Optimus", work_space_json)

        if updated_workflow.status_code != 200:
            print("ERROR :" + updated_workflow.content)
            sys.exit(1)
        else:
            print("updated successfully")

    elif args.cmd == 'submit_workflow':
        # Launching the Updated Monitor Submission Workflow
        response = fapi.get_entities_with_type(billing_project, args.workspace_name)
        entities = response.json()

        for ent in entities:
            ent_name = ent['name']
            ent_type = ent['entityType']
            ent_attrs = ent['attributes']

        submisson_response = fapi.create_submission(billing_project, args.workspace_name,
                                                    args.workflow_repo, args.workflow_name,
                                                    entity=args.entity_id, etype="participant_lane_set",
                                                    expression=None, use_callcache=True)

        if submisson_response.status_code != 201:
            print(submisson_response.content)
            sys.exit(1)
        else:
            print("Successfully Created Submisson")
            with open('response.txt', 'w') as fout:
                # json.dump(submisson_response.json(), fout)
                fout.write(submisson_response.json()['submissionId'] + '\n')
        # r = create_workspace_config("broadgdac", args.workspace_name, body):
        # print(r.content.decode())
    elif args.cmd == 'get_status':
        res = fapi.get_submission(billing_project, args.workspace_name, args.submission_id)
        print(res.content.decode())


if __name__ == "__main__":
    main()
