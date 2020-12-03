#!/usr/bin/env python3

import argparse
import loompy
import numpy as np
import tempfile
import os


def main():
    description = """Combine library level loom files into a single project level loom and add global metadata. 
    Cell barcodes from separate libraries are suffixed with a number to avoid collisions."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-loom-files',
                        dest='input_loom_files',
                        nargs="+",
                        required=False,
                        help="Paths to input loom files")
    parser.add_argument('--library',
                        dest='library',
                        nargs="+",
                        required=False,
                        help="Library metadata string")
    parser.add_argument('--species',
                        dest='species',
                        nargs="+",
                        required=False,
                        help="Species metadata string")
    parser.add_argument('--organ',
                        dest='organ',
                        nargs="+",
                        required=False,
                        help="Organ metadata string")
    parser.add_argument('--project-id',
                        dest='project_id',
                        required=False,
                        help="Project ID")
    parser.add_argument('--project-name',
                        dest='project_name',
                        required=False,
                        help="Project Name")
    parser.add_argument('--output-loom-file',
                        dest='output_loom_file',
                        required=False,
                        help="Path to output loom file")

    parser.add_argument('--run-tests',
                        dest='run_tests',
                        required=False,
                        action = 'store_true',
                        help="Runs tests")

    args = parser.parse_args()

    if args.run_tests:
       run_tests()
    else:
      if args.project_id:
          combine_loom_files(loom_file_list = args.input_loom_files,
                             library = ", ".join(set(args.library)),
                             species = ", ".join(set(args.species)),
                             organ = ", ".join(set(args.organ)),
                             project_id = args.project_id,
                             project_name = args.project_name,
                             output_loom_file = args.output_loom_file
                            )

def run_tests():
    # merging three files
    with tempfile.TemporaryDirectory() as tmpdirname:
        print('created temporary directory', tmpdirname)

        input_loom_files = [ 'a.loom', 'b.loom', 'c.loom' ] 
        input_libraries = [ 'V2', 'V2', 'V2']
        input_species = [ 'human', 'human', 'human']
        input_organs = [ 'liver', 'lung', 'brain']
        input_project_id = "project_id"
        input_project_name = "project_name"
        output_loom_file = "output_loom_file.loom"
        input_loom_ncols = [ 10, 10, 10 ]
        input_ids = [ 'CODE1', 'CODE2', 'CODE3' ]
        input_names = [ 'NAME1', 'NAME2', 'NAME3' ]

        for input_id, input_name, filename, ncol in \
            zip(input_ids, input_names, input_loom_files, input_loom_ncols): 

            matrix = np.arange(ncol*10).reshape(10, ncol)
            gene_names =  [ 'gene_' + str(ord('a') + i) for i in range(10) ]

            row_attrs = { "row_attribute_1": np.arange(10), 
                          "row_attribute_2": np.arange(10),
                          "gene_names": gene_names
                        }

            cell_names =  [ 'cell_' +  str(ord('a') + i) for i in range(10) ]
            col_attrs = { "column_attribute_1": np.arange(ncol) ,
                          "cell_names" : cell_names,
                          "n_molecules": np.arange(95, 95 + ncol)
                        }

            loompy.create(filename, matrix, row_attrs, col_attrs)  

            ds = loompy.connect(filename)
            ds.attrs["expression_data_type"] = 'exonic'
            ds.attrs["optimus_output_schema_version"] = '1.0.0'
            ds.attrs["pipeline_version"] = 'Optimus_v4.1.7'
            ds.attrs["input_id_metadata_field"] = 'sequencing_process.provenance.document_id'
            ds.attrs["input_name_metadata_field"] = 'sequencing_input.biomaterial_core.biomaterial_id'
            ds.attrs["input_id"] = input_id
            ds.attrs["input_name"] = input_name

            ds.close()

            

        combine_loom_files(loom_file_list = input_loom_files,
                             library = ", ".join(set(input_libraries)),
                             species = ", ".join(set(input_species)),
                             organ = ", ".join(set(input_organs)),
                             project_id = input_project_id,
                             project_name = input_project_name,
                             output_loom_file = output_loom_file
                           )

        with loompy.connect(output_loom_file) as dsin: 
            expected_cell_names =  ['cell_102-0', 'cell_103-0', 'cell_104-0', 'cell_105-0', 'cell_106-0', 
                               'cell_102-1', 'cell_103-1', 'cell_104-1', 'cell_105-1', 'cell_106-1', 
                               'cell_102-2', 'cell_103-2', 'cell_104-2', 'cell_105-2', 'cell_106-2']

            expected_gene_names = ['gene_100', 'gene_101', 'gene_102', 'gene_103', 
                                   'gene_104', 'gene_105', 'gene_106', 'gene_97', 
                                  'gene_98', 'gene_99']

            # test the expected number of barcodes after filtering
            assert(sorted(list(dsin.ca['cell_names'])) == sorted(expected_cell_names))

            # check the number of gene-names are as expected
            assert(sorted(list(dsin.ra['gene_names'])) == sorted(expected_gene_names))
        
            # the expected input_id in the merged loom must be concatenated
            assert(dsin.attrs["input_id"]  == ", ".join(input_ids))

            # the expected input_names in the merged loom must be concatenated
            assert(dsin.attrs["input_name"] == ", ".join(input_names))

            # the expression data type should be the same
            assert(dsin.attrs["expression_data_type"] == 'exonic')



def combine_loom_files(loom_file_list, library, species, organ, project_id, project_name, output_loom_file):
    expression_data_type_list = []
    optimus_output_schema_version_list = []
    pipeline_versions_list = []
    input_id_metadata_field_list = []
    input_name_metadata_field_list = []
    input_id_list = []
    input_name_list = []

    with loompy.new(output_loom_file) as dsout:
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
    ds = loompy.connect(output_loom_file)

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

