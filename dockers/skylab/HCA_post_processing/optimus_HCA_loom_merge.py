#!/usr/bin/env python3

import argparse
import loompy


def main():
    description = """Add metadata into a Loom file as column attributes"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-loom-files',
                        dest='input_loom_files',
                        nargs="+",
                        required=True,
                        help="Paths to input loom files")
    parser.add_argument('--library',
                        dest='library',
                        nargs="+",
                        required=True,
                        help="Library metadata string")
    parser.add_argument('--species',
                        dest='species',
                        nargs="+",
                        required=True,
                        help="Species metadata string")
    parser.add_argument('--stage',
                        dest='stage',
                        nargs="+",
                        required=True,
                        help="Stage metadata string")
    parser.add_argument('--organ',
                        dest='organ',
                        nargs="+",
                        required=True,
                        help="Organ metadata string")
    parser.add_argument('--output-loom-file',
                        dest='output_loom_file',
                        required=True,
                        help="Path to output loom file")
    args = parser.parse_args()

    loom_file_list = args.input_loom_files
    library = " ".join(args.library)
    species = " ".join(args.species)
    stage = " ".join(args.stage)
    organ = " ".join(args.organ)

    attr_dict = {"library_costruction_method": library,
                 "species": species,
                 "devleopmental_stage": stage,
                 "organ": organ
                 }

    expression_data_type_list = []
    optimus_output_schema_version_list = []
    pipeline_versions_list = []

    # To avoid collisions, the barcoodes for each loom file will have a numerical extension added
    for i in range(len(loom_file_list)):
        loom_file = loom_file_list[i]
        ds = loompy.connect(loom_file)

        # add the file index as an extension to the barcode to ensure ther are no collisions
        ds.ca['cell_names'] = ds.ca['cell_names'] + "-" + str(i)

        # add input_id and input_name as column attributes
        num_rows, num_cols = ds.shape
        input_id = ds.attrs['input_id']
        ds.ca.input_id = [input_id for x in range(num_cols)]
        try:
            input_name = ds.attrs['input_name']
            ds.ca.input_name = [input_name for x in range(num_cols)]
        except AttributeError:
            pass

        # add global attributes for this file to the running list of global attributes
        expression_data_type_list.append(ds.attrs['expression_data_type'])
        optimus_output_schema_version_list.append(ds.attrs['optimus_output_schema_version'])
        pipeline_versions_list.append(ds.attrs['pipeline_version'])

        ds.close()

    attr_dict['expression_data_type'] = ", ".join(set(expression_data_type_list))
    attr_dict['optimus_output_schema_version'] = ", ".join(set(optimus_output_schema_version_list))
    attr_dict['pipeline_version'] = ", ".join(set(pipeline_versions_list))

    # comobine the loom files
    loompy.combine(loom_file_list, output_file=args.output_loom_file, file_attrs=attr_dict)

    # alter the global attributes of the combired loom file
    ds = loompy.connect(args.output_loom_file)

    # delete the input_id and input_name (which are now column attributes)
    del ds.attrs.input_id
    try:
        del ds.attrs.input_name
    except KeyError:
        pass

    ds.close()


if __name__ == '__main__':
    main()

