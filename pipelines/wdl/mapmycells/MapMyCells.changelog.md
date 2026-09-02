# 0.1.0
2026-09-02 (Date of Last Commit)

* Added the MapMyCells pipeline, wrapping the Allen Institute `cell_type_mapper` library for cell type label transfer against precomputed reference atlases (Human_MTG, Mouse_WMB, or a Custom taxonomy)
* Moved the Human_MTG / Mouse_WMB precomputed-stats and marker-gene assets into the docker image (built in warp-tools) instead of referencing `gs://` paths that were never actually populated; the task now selects the in-container path per `reference_atlas` instead of computing a workflow-level `File` from a bucket that didn't have the file
* Added an `ErrorWithMessage` check that fails fast when `reference_atlas` is `"Custom"` but `custom_precomputed_stats` is not supplied, instead of failing later with an opaque `select_first` error
* Repinned the `docker` default to the image tag actually published for this branch (`add-mapmycells`); the previous default tag was never built
* Repointed `gene_mapping_db`'s default (applies to every `reference_atlas`, not just Human_MTG/Mouse_WMB) at a real, now-populated copy of Allen Institute's pre-built gene-symbol-to-Ensembl-ID db at `gs://broad-gotc-test-storage/mapmycells/mmc_gene_mapper.2025-08-04.db` (~15GB); the previous path pointed at a file that never existed. Kept this asset GCS-hosted rather than baked into the docker image -- warp-tools has no precedent for baking reference data this large into an image, and doing so would balloon every task's docker pull by ~15GB regardless of whether that run uses gene mapping. Bumped `disk_size` default 100->150 to cover its localization
* Overridable per-run via the workflow-level `custom_gene_mapping_db` input
