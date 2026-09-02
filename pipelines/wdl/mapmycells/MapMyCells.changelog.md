# 0.1.0
2026-09-02 (Date of Last Commit)

* Added the MapMyCells pipeline, wrapping the Allen Institute `cell_type_mapper` library for cell type label transfer against precomputed reference atlases (Human_MTG, Mouse_WMB, or a Custom taxonomy)
* Moved the Human_MTG / Mouse_WMB precomputed-stats and marker-gene assets into the docker image (built in warp-tools) instead of referencing `gs://` paths that were never actually populated; the task now selects the in-container path per `reference_atlas` instead of computing a workflow-level `File` from a bucket that didn't have the file
* Added an `ErrorWithMessage` check that fails fast when `reference_atlas` is `"Custom"` but `custom_precomputed_stats` is not supplied, instead of failing later with an opaque `select_first` error
* Repinned the `docker` default to the image tag actually published for this branch (`add-mapmycells`); the previous default tag was never built
* Dropped the `gene_mapping_db` default for Human_MTG/Mouse_WMB (it pointed at a `gs://` file that didn't exist); gene-symbol-to-Ensembl mapping is now opt-in only via `custom_gene_mapping_db`, for every atlas, pending a decision on baking in `mmc_gene_mapper`'s db (its build step downloads 15GB of NCBI taxonomy data)
