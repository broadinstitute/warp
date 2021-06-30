# Consortia Data Processing

## Human Cell Atlas Data Coordination Platform Matrix Processing
Optimus supports data processing for the Human Cell Atlas (HCA) Data Coordination Platform (DCP). Learn more about the DCP at the [HCA Data Portal](https://data.humancellatlas.org/)).

All DCP Projects processed with Optimus have matrices containing the standard metrics and counts detailed in the [Optimus Count Matrix Overview](./Loom_schema.md), but also have additional post-processing to incorporate DCP-curated metadata.

**To reduce file size, DCP Loom matrices are in sparse format and are minimally filtered so that only cells with 100 molecules or more are retained.**

:::warning
This section details matrices produced for the Human Cell Atlas (HCA) [Data Coordination Platform (DCP)2.0](https://data.humancellatlas.org/), which includes matrices processed with Optimus v4.1.7 and later. The DCP is currently reprocessing data generated with earlier Optimus versions and will deprecate previous matrices once reprocessing is complete.
:::

DCP matrices contain DCP-curated metadata in the Loom global attributes (see table below) which may be useful when exploring the data and linking it back to the Project metadata. Read more about each metadata field in the DCP [Metadata Dictionary](https://data.humancellatlas.org/metadata).

| Metadata Attribute Name in Count Matrix | Metadata Description |
| --- | --- |
| `donor_organism.genus_species` | species information; human or mouse |
| `library_preparation_protocol.library_construction_approach` | technology used for library preparation, i.e 10x or SS2 |
| `specimen_from_organism.organ` | organ |
| `project.project_core.project_name` | project name |
| `project.provenance.document_id` | project id |
| `input_id` | metadata values for  `sequencing_process.provenance.document_id`; unique ID to demarcate the library prep |
| `input_name` | metadata values for `sequencing_input.biomaterial_core.biomaterial_id`; unique ID for the biomaterial |
| `input_id_metadata_field` | string describing the DCP-curated metadata field used for input_id: `sequencing_process.provenance.document_id` |
| `input_name_metadata_field` | string describing the DCP-curated metadata field used for input_name: `sequencing_input.biomaterial_core.biomaterial_id` |

To create the DCP project matrices, Loom outputs from individual 10x library preparations, each with their own `input_id`, are combined into a single Loom file.


### Making matrix cell barcodes unique for DCP project matrices
Since DCP project matrices often contain combined data from multiple library preparations, the project matrix cell bardcodes are modified so that they are unique for each library preparation, allowing the barcodes to be used by downstream community tools like Cumulus and Seurat.

In the generic Optimus matrices, cell barcodes are listed in both the the `CellID` and the `cell_names` columns.

For DCP projects, however, **the cell barcodes in the  `cell_names` column are modified** so that each cell barcode belonging to an individual library preparation is unique. This is done by adding a numerical suffix to the barcodes that corresponds to `input_id` for the library preparation from which the cell barcodes came. 

The `input_ids` are listed in the matrix global attributes. The order of the input_ids serves as index for the cell barcode suffix.

For example, let's look at the global attribute `input_id` for a DCP project matrix:

```python
>>> ds.attrs.input_id
'166c1b1a-ad9c-4476-a4ec-8b52eb5032c7, 22b7da3d-a301-433e-99e1-e67266c1ee8b, 337a48c5-e363-45aa-886f-ccd4425edc2b, 40630e8b-c3a3-4813-b1e4-b156637c5cc3, 58d703d1-d366-42d0-af44-a3bb836838a5, 70c8d647-7984-4d03-912a-f2437aa1ba4f, 7c86cf30-4284-4a0d-817f-6047560c05c3, 8ef7aca4-be00-4c03-8576-1b2eff4ce7af, ae0cfa6e-e7cb-4a88-9f89-1c44abaa2291, cbd23025-b1bf-4e9e-a297-ddab4a217b76, df049da4-3d20-4da7-a1d7-7d6e8f7740ff, e17bf5ea-788b-4756-a008-a07aec091e10'
>>> 
```
Each of these UUIDs represents one library preparation. This matrix contains data from 12 library preparations total. 

Now let's look at the `cell_names` column attribute which contains the unique cell barcodes:

```python
>>> ds.ca.cell_names
array(['GGACAAGAGTGCGTGA-0', 'GATCGATCACCAGGTC-0', 'AGCGGTCAGGGCTTGA-0',
       ..., 'GTACGTAAGCTATGCT-11', 'CAGAATCTCTGAGTGT-11',
       'AACACGTAGTGTTTGC-11'], dtype=object)
>>> 
```

The suffix appended to the barcodes in the `cell_names` column is the index for the `input_id` UUID to which the cell barcodes belong. 

For example, cell barcodes with a "-0" suffix belong to the library preparation represented by the first UUID, `166c1b1a-ad9c-4476-a4ec-8b52eb5032c7`, whereas cell barcodes with a "-11" suffix represent the 10th UUID, `cbd23025-b1bf-4e9e-a297-ddab4a217b76`.


### Mapping DCP project matrix data to the metadata manifest

While the project matrices contain some project metadata (listed in the table above), there is additionally useful metadata in the project metadata manifest, a TSV file containing all of a project's metadata, including donor and disease state information.

The project matrix `input_id` can be used to map the data back to the DCP metadata manifest. Each manifest contains a column called `sequencing_process.provenance.document_id`. The values in this column will match the matrix `input_id`.

Read more about the metadata manifest in the DCP [Exploring Projects guide](https://data.humancellatlas.org/guides).

:::tip Explore HCA Project matrices in Terra
HCA matrices produced with Optimus are compatible with multiple downstream community analysis tools. For a tutorial on using the Optimus matrix with [Seurat](https://satijalab.org/seurat/index.html), [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Cumulus](https://cumulus.readthedocs.io/en/latest/index.html), or [Pegasus](https://pegasus.readthedocs.io/en/stable/#), see the public [Intro-to-HCA-data-on-Terra workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/Intro-to-HCA-data-on-Terra) (login required) and its accompanying [step-by-step guide](https://support.terra.bio/hc/en-us/articles/360060041772).