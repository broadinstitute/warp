#!/usr/bin/env python3
from optimus_HCA_loom_merge import combine_loom_files
import loompy
import numpy as np
import os
import re
import pytest

_test_data_dir = os.path.join(os.path.split(__file__)[0], "testdata/")


def test_combine_loom_files(tmpdir):
    import os
    print('Test to merge three loom files')
    # original file 
    looma = os.path.join(_test_data_dir, "a.loom")
    loomb = os.path.join(_test_data_dir, "b.loom")
    loomc = os.path.join(_test_data_dir, "c.loom")
    loomd = os.path.join(_test_data_dir, "d.loom")
    
    input_loom_files = [looma, loomb, loomc]
    input_libraries = ['V2', 'V2', 'V2']
    input_species = ['human', 'human', 'human']
    input_organs = ['liver', 'lung', 'brain']
    input_project_id = "project_id"
    input_project_name = "project_name"
    output_loom_file = "output_loom_file.loom"

    combine_loom_files(loom_file_list=input_loom_files,
                       library=", ".join(set(input_libraries)),
                       species=", ".join(set(input_species)),
                       organ=", ".join(set(input_organs)),
                       project_id=input_project_id,
                       project_name=input_project_name,
                       output_loom_file=output_loom_file)

    # open the loom files
    dsina = loompy.connect(looma)
    dsinb = loompy.connect(loomb)
    dsinc = loompy.connect(loomc)
    dsin = loompy.connect(output_loom_file)

    # the expected input_id in the merged loom must be concatenated
    input_ids = ['58a18a4c-5423-4c59-9b3c-50b7f30b1ca5',
                 'c763f679-e13d-4f81-844f-c2c80fc90f46',
                 'c76d90b8-c190-4c58-b9bc-b31f586ec7f2']
    assert(dsin.attrs["input_id"] == ", ".join(input_ids))

    # the expected input_names in the merged loom must be concatenated
    input_names = ['PP012_suspension', 'PP003_suspension', 'PP004_suspension']
    assert(dsin.attrs["input_name"] == ", ".join(input_names))

    # the expression data type should be the same
    assert(dsin.attrs["expression_data_type"] == 'exonic')

    assert(dsin.attrs["library_preparation_protocol.library_construction_approach"]
           == ", ".join(set(input_libraries)))
    assert(dsin.attrs["donor_organism.genus_species"] == ", ".join(set(input_species)))
    assert(dsin.attrs["specimen_from_organism.organ"] == ", ".join(set(input_organs)))
    assert(dsin.attrs["project.provenance.document_id"] == input_project_id)
    assert(dsin.attrs["project.project_core.project_name"] == input_project_name)

    # the number of cells should be 30
    assert(dsin.shape == (58347, 30))

    # the list of gene names must be the same in all the loom files
    assert(all(dsin.ra['gene_names'] == dsina.ra['gene_names']))
    assert(all(dsin.ra['gene_names'] == dsinb.ra['gene_names']))
    assert(all(dsin.ra['gene_names'] == dsinc.ra['gene_names']))

    # the output loom file must have the attribute keys in the individual looms
    attrs_set = set(list(dsin.attrs.keys()))
    attrs_seta = set(list(dsina.attrs.keys()))
    attrs_setb = set(list(dsinb.attrs.keys()))
    attrs_setc = set(list(dsinc.attrs.keys()))
    assert(attrs_seta.difference(attrs_set) == set())
    assert(attrs_setb.difference(attrs_set) == set())
    assert(attrs_setc.difference(attrs_set) == set())

    # test the counts for the cells in a.loom
    sample_no_regex = re.compile(r'-(\d+)$')
    for i, barcode in enumerate(dsin.ca['cell_names']):
        res = sample_no_regex.search(barcode)
        if res:
            sample_no = res.group(1)
        if sample_no == "0":
            original_barcode = re.sub(r'-0$', '', barcode)
            dsin_o = dsina
        elif sample_no == "1":
            original_barcode = re.sub(r'-1$', '', barcode)
            dsin_o = dsinb
        elif sample_no == "2":
            original_barcode = re.sub(r'-2$', '', barcode)
            dsin_o = dsinc
        else:
            # no barcode should fail all the preceeding predicates
            assert False

        mismatched_indices = np.where(dsin[:, np.where(dsin.ca['cell_names'] == barcode)[0]] !=
                                      dsin_o[:, np.where(dsin_o.ca['cell_names'] == original_barcode)[0]])

        assert(len(mismatched_indices[0]) == 0)

        if len(mismatched_indices[0]):
            index = np.where(dsin.ca['cell_names'] == barcode)[0]
            index_o = np.where(dsin_o.ca['cell_names'] == original_barcode)[0]
            print()
            vals = [dsin[i, index][0] for i in mismatched_indices[0]]
            vals_o = [dsin_o[i, index_o][0] for i in mismatched_indices[0]]
            print(index, index_o, dsin_o.ca['cell_names'][index_o][0], dsin.ca['cell_names'][index][0])
            print(dsin.ra['gene_names'][mismatched_indices[0]])
            print(dsin_o.ra['gene_names'][mismatched_indices[0]])
            print(vals)
            print(vals_o)

        with pytest.raises(AssertionError):
            combine_loom_files(loom_file_list=[looma, loomd],
                               library="10X 3' v2 sequencing",
                               species="Homo sapiens",
                               organ="brain",
                               project_id=input_project_id,
                               project_name=input_project_name,
                               output_loom_file="invalid_matrix.loom")

    dsin.close()
    dsina.close()
    dsinb.close()
    dsinc.close()
