import argparse
import loompy
import numpy as np
import scipy as sc


def convert_to_sparse_matrix(count_matrix, nrows, ncols):
    # read loom file expression counts and add it to the expression_counts dataset
    exp_counts = np.load(count_matrix)
    # now convert it back to a csr_matrix object
    csr_exp_counts = sc.sparse.csr_matrix(
        (exp_counts["data"], exp_counts["indices"], exp_counts["indptr"]),
        shape=exp_counts["shape"],
    )

    expr_sp = sc.sparse.coo_matrix((nrows, ncols), np.float32)

    xcoord = []
    ycoord = []
    value = []

    chunk_row_size = 10000
    chunk_col_size = 10000

    for i in range(0, nrows, chunk_row_size):
        for j in range(0, ncols, chunk_col_size):
            p = chunk_row_size
            if i + chunk_row_size > nrows:
                p = nrows - 1
            q = chunk_col_size
            if j + chunk_col_size > ncols:
                q = ncols - j
            expr_sub_row_coo = sc.sparse.coo_matrix(csr_exp_counts[i:i + p, j:j + q].toarray())
            for k in range(0, expr_sub_row_coo.data.shape[0]):
                xcoord.append(expr_sub_row_coo.row[k] + i)
                ycoord.append(expr_sub_row_coo.col[k] + j)
                value.append(expr_sub_row_coo.data[k])

    xcoord = np.asarray(xcoord)
    ycoord = np.asarray(ycoord)
    value = np.asarray(value)

    expr_sp_t = sc.sparse.coo_matrix((value, (ycoord, xcoord)), shape=(expr_sp.shape[1], expr_sp.shape[0]))

    del xcoord
    del ycoord
    del value

    return expr_sp_t


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

    combine_loom_files(loom_file_list = args.input_loom_files,
                             library = ", ".join(set(args.library)),
                             species = ", ".join(set(args.species)),
                             organ = ", ".join(set(args.organ)),
                             project_id = args.project_id,
                             project_name = args.project_name,
                             output_loom_file = args.output_loom_file
                            )

def combine_loom_files(loom_file_list, library, species, organ, project_id, project_name, output_loom_file):
    expression_data_type_list = []
    optimus_output_schema_version_list = []
    pipeline_versions_list = []
    input_id_metadata_field_list = []
    input_name_metadata_field_list = []
    input_id_list = []
    input_name_list = []

    with loompy.new("intermediate.loom") as dsout:
        for i in range(len(loom_file_list)):
            loom_file = loom_file_list[i]
            with loompy.connect(loom_file) as ds:

                # add global attributes for this file to the running list of global attributes
                expression_data_type_list.append(ds.attrs["expression_data_type"])
                optimus_output_schema_version_list.append(ds.attrs["optimus_output_schema_version"])
                pipeline_versions_list.append(ds.attrs["pipeline_version"])
                input_id_metadata_field_list.append(ds.attrs["input_id_metadata_field"])
                input_name_metadata_field_list.append(ds.attrs["input_name_metadata_field"])
                input_id_list.append(ds.attrs["input_id"])
                input_name_list.append(ds.attrs["input_name"])

                # check that the ordering is the same for the matrices being combined
                if dsout.shape[0] != 0:
                    assert(np.array_equal(dsout.ra["ensembl_ids"], ds.ra["ensembl_ids"]))

                # filter out cells with low counts n_molecules > 1
                UMIs = ds.ca['n_molecules']
                cells = np.where(UMIs >= 100)[0]
                for (ix, selection, view) in ds.scan(items=cells, axis=1):
                    view.ca['cell_names'] = view.ca['cell_names'] + "-" + str(i)
                    dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

    # add global attributes for this file to the running list of global attributes
    ds = loompy.connect("intermediate.loom")

    # Grab col attributes and row attributes from the merged matrix
    row_attrs = ds.ra[:]
    col_attrs = ds.ca[:]

    # Grab count matrix from merged matirx
    count_matrix = ds[:, :]
    nrows, ncols = ds.shape
    ds.close()

    # Convert the count matirx to a sparse matrix
    sparse_matrix = convert_to_sparse_matrix(count_matrix, nrows, ncols)

    # Write out a new loom file with the sparse matrix
    loompy.create(output_loom_file, sparse_matrix, row_attrs, col_attrs)

    # add the global atributes to the loom file
    ds = loompy.connect(output_loom_file)

    ds.attrs["library_preparation_protocol.library_construction_approach"] = library
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

