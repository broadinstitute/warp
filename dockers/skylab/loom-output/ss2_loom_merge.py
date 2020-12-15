#!/usr/bin/env python3

import argparse
import loompy


def main():
    description = """Merge the outputs of multiple SS2 pipeline runs into a single Loom file"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-loom-files',
                        dest='input_loom_files',
                        nargs="+",
                        required=True,
                        help="Path to input loom directory in DirectoryStore format")
    parser.add_argument('--output-loom-file',
                        dest='output_loom_file',
                        required=True,
                        help="Path to output loom file")
    parser.add_argument('--batch_id',
                        dest='batch_id',
                        required=True,
                        help="Batch id for output loom")
    parser.add_argument('--batch_name',
                        dest='batch_name',
                        help='User provided plate id for output loom')
    parser.add_argument('--project_id',
                        dest='project_id',
                        help='User provided plate id for output loom')
    parser.add_argument('--project_name',
                        dest='project_name',
                        help='User provided plate name for output loom')
    parser.add_argument('--library',
                        dest='library',
                        help='User provided library for output loom')
    parser.add_argument('--species',
                        dest='species',
                        help='User provided species for output loom')
    parser.add_argument('--organ',
                        dest='organ',
                        help='User provided organ for output loom')
    parser.add_argument('--pipeline_version',
                        dest='pipeline_version',
                        required=True,
                        help='Multisample SS2 version')
    args = parser.parse_args()

    # The list of Loom files that we need to merge

    loom_file_list = args.input_loom_files

    attrDict = dict()
    attrDict['batch_id'] = args.batch_id
    attrDict['pipeline_version'] = args.pipeline_version
    if args.batch_name is not None:
        attrDict['batch_name'] = args.batch_name

    if args.library is not None:
        attrDict['library_preparation_protocol.library_construction_approach'] = args.library

    if args.species is not None:
        attrDict['donor_organism.genus_species'] = args.species

    if args.organ is not None:
        attrDict['specimen_from_organism.organ'] = args.organ

    if args.project_id is not None:
        attrDict['project.provenance.document_id'] = args.project_id

    if args.project_name is not None:
        attrDict['project.project_core.project_short_name'] = args.project_name

    loompy.combine(loom_file_list,output_file=args.output_loom_file,file_attrs = attrDict)

    #TODO: check global attributes and make sure that they are correct and make sense (i.e. no input_id)


if __name__ == '__main__':
    main()
