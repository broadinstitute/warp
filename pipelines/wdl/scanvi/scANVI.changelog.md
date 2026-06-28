# 1.1.0
2026-06-27 (Date of Last Commit)

* Made the ATAC h5ad input optional. When no ATAC h5ad is provided (direct-file mode) or no ATAC object is found in the input bucket (bucket mode), the pipeline auto-detects GEX-only mode and trains/annotates from the reference atlas using gene expression and reference data alone, without using ATAC.
* Added run_gex_only_model to the scvi-scanvi container to train SCVI/SCANVI on GEX + reference only
* PreprocessFilter skips ATAC loading, barcode reindexing, shared-barcode subsetting, and gene-activity conversion in GEX-only mode
* Made preprocessed_atac_activity_h5ad and atac_annotated_h5ad optional outputs (null in GEX-only mode)
* Added a GEX-only example inputs file (example_inputs/scANVI.gex_only.json)
* Added an optional max_epochs input that caps SCVI/SCANVI training epochs in both multiome and GEX-only modes (defaults to the container value of 500 when unset)

# 1.0.0
2026-04-17 (Date of Last Commit)

* Initial release of the ScviScanvi pipeline for cell type label transfer on Multiome data
* Integrated SCVI and SCANVI models using the scvi-scanvi docker image (1.0.0-1.2-1756234975)
* Test GPU support for accelerated model training
* Outputs SCANVI predictions, annotated ATAC, and annotated GEX h5ad files
* Added PreprocessFilter task to handle all h5ad preprocessing and filtering as a separate CPU-only step before GPU model training
* Moved column patching (star_IsCell, gex_barcodes), GEX cell/gene filtering, ATAC barcode reindexing, shared barcode subsetting, batch/modality tagging, and gene activity matrix conversion out of the monolithic MultiomeLabelTransfer task
* Simplified MultiomeLabelTransfer task to accept preprocessed h5ad files directly; removed inline Python patching and bucket resolution logic
* Hardcode GPU configuration in runtime block for MultiomeLabelTransfer; no GPU needed for PreprocessFilter
* added test_cuda.wdl to the tasks directory

