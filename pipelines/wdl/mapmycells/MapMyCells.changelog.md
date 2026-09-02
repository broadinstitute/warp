# 0.1.0
2026-09-02 (Date of Last Commit)

* Added the MapMyCells pipeline, wrapping the Allen Institute `cell_type_mapper` library for cell type label transfer against precomputed reference atlases (Human_MTG, Mouse_WMB, or a Custom taxonomy)
* Fixed the `query_markers` selection logic so the Custom/Human_MTG fallback is a single explicit condition instead of a dead branch that duplicated the same value on both sides of the conditional
* Added an `ErrorWithMessage` check that fails fast when `reference_atlas` is `"Custom"` but `custom_precomputed_stats` is not supplied, instead of failing later with an opaque `select_first` error
* Repinned the `docker` default to the image tag actually published for this branch (`add-mapmycells`); the previous default tag was never built
