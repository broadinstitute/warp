# MultiSampleSmartSeq2_v2.1.0

2020-08-10 (Date of Last Commit)

* Added batch_name as an optional input for user provided biomaterial id.
* Passed pipeline_version to output loom file  
* Added input_id_metadata_field and input_name_metadata_field as optional input

# MultiSampleSmartSeq2_v2.0.1

2020-07-20 (Date of Last Commit)

* Changed the imports to relative imports to support Dockstore->Terra release

# MultiSampleSmartSeq2_v2.0.0

2020-06-04 (Date of Last Commit)

* Removed zarr output and made loom output as default

* Loom file attribute names have changed: CellID: cell_names, Gene: gene_names and Accession: ensembl_ids

* Loom file name has changed from out.loom to "plateid".loom

* Added the expected counts in addition to TPMs in the loom matrix. 

# MultiSampleSmartSeq2_v1.1.0

2020-05-07 (Date of Last Commit)

* Added estimated count matrix to zarr output

# MultiSampleSmartSeq2_v1.0.0

2019-12-16 (Date of Last Commit)

* This is the first release of the Smart-seq2 Multi Sample workflow
