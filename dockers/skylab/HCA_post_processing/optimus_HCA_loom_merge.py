#!/usr/bin/env python3

import argparse
import loompy
import numpy as np


def main():
    description = """Combine library level loom files into a single project level loom and add global metadata. 
    Cell barcodes from separate libraries are suffixed with a number to avoid collisions."""
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
    parser.add_argument('--organ',
                        dest='organ',
                        nargs="+",
                        required=True,
                        help="Organ metadata string")
    parser.add_argument('--project-id',
                        dest='project_id',
                        required=True,
                        help="Project ID")
    parser.add_argument('--project-name',
                        dest='project_name',
                        required=True,
                        help="Project Name")
    parser.add_argument('--output-loom-file',
                        dest='output_loom_file',
                        required=True,
                        help="Path to output loom file")
    args = parser.parse_args()

    loom_file_list = args.input_loom_files
    library = ", ".join(set(args.library))
    species = ", ".join(set(args.species))
    organ = ", ".join(set(args.organ))
    project_id = args.project_id
    project_name = args.project_name

    expression_data_type_list = []
    optimus_output_schema_version_list = []
    pipeline_versions_list = []
    input_id_metadata_field_list = []
    input_name_metadata_field_list = []
    input_id_list = []
    input_name_list = []

    with loompy.new(args.output_loom_file) as dsout:
        for i in range(len(loom_file_list)):
            loom_file = loom_file_list[i]
            with loompy.connect(loom_file) as ds:
                ds.ca['cell_names'] = ds.ca['cell_names'] + "-" + str(i)
                # add input_id and input_name as column attributes
                # num_rows, num_cols = ds.shape
                # input_id = ds.attrs['input_id']
                # ds.ca.input_id = [input_id for x in range(num_cols)]
                # try:
                    # input_name = ds.attrs['input_name']
                    # ds.ca.input_name = [input_name for x in range(num_cols)]
                # except AttributeError:
                    # pass

                # add global attributes for this file to the running list of global attributes
                expression_data_type_list.append(ds.attrs["expression_data_type"])
                optimus_output_schema_version_list.append(ds.attrs["optimus_output_schema_version"])
                pipeline_versions_list.append(ds.attrs["pipeline_version"])
                input_id_metadata_field_list.append(ds.attrs["input_id_metadata_field"])
                input_name_metadata_field_list.append(ds.attrs["input_name_metadata_field"])
                input_id_list.append(ds.attrs["input_id"])
                input_name_list.append(ds.attrs["input_name"])

                # filter out cells with low counts n_molecules > 1
                UMIs = ds.ca['n_molecules']
                cells = np.where(UMIs >= 100)[0]
                for (ix, selection, view) in ds.scan(items=cells, axis=1, key="gene_names"):
                    dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

    # add global attributes for this file to the running list of global attributes
    ds = loompy.connect(args.output_loom_file)

    ds.attrs["library_preparation_protocol.library_construction_method"] = library
    ds.attrs["donor_organism.genus_species"] = species
    ds.attrs["specimen_from_organism.organ"] = organ
    ds.attrs["project.provenance.document_id"] = project_id
    ds.attrs["project.project_core.project_name"] = project_name
    ds.attrs["expression_data_type"] = ", ".join(set(expression_data_type_list))
    ds.attrs["optimus_output_schema_version"] = ", ".join(set(optimus_output_schema_version_list))
    ds.attrs["pipeline_version"] = ", ".join(set(pipeline_versions_list))
    ds.attrs["input_id_metadata_field"] = ", ".join(set(input_id_metadata_field_list))
    ds.attrs["input_name_metadata_field"] = ", ".join(set(input_name_metadata_field_list))
    ds.attrs["input_id"] = ", ".join(input_id_list)
    ds.attrs["input_name"] = ", ".join(input_name_list)

    ds.close()


if __name__ == '__main__':
    main()

